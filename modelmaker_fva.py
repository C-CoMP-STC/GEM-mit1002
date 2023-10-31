# modelmaker_fva.py
'''
Slight modification of Cobrapy's implementation of FVA that uses pebble instead
of multiprocessing when running in parallel because for whatever reason pebble
tolerates reactions with infinite bounds and Cobrapy's multiprocessing doesn't
'''

import pandas as pd
import numpy as np
from optlang.interface import Objective
from optlang.symbolics import Zero
from warnings import warn
from cobra.flux_analysis.parsimonious import add_pfba
import itertools as it
from pebble import ProcessPool, ProcessExpired
import time
from concurrent.futures import TimeoutError
from optlang.interface import OPTIMAL, UNBOUNDED

def fva(
    model, reaction_list = None, fraction_of_optimum = 1.0, pfba_factor = None,
    processes = 1, verbose = 0, update_freq = 1000
):
    """Determine the minimum and maximum flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default None).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum (default 1.0).
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the `pfba_factor`
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds
        (default None).
    processes : int, optional
        The number of parallel processes to run. 1 by default

    Returns
    -------
    pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes
    -----
    This implements the fast version as described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    sub-optimal.

    References
    ----------
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235
    """
    # use all reactions in model if none were given
    if reaction_list is None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        if any(not isinstance(r, str) for r in reaction_list):
            msg = 'All elements of reaction_list must be IDs (i.e. strings), '
            msg += 'not cobra.Reaction objects or anything else.'
            raise ValueError(msg)
        else:
            reaction_ids = reaction_list
    # make sure we're not using more processes than reactions
    num_reactions = len(reaction_ids)
    processes = min(processes, num_reactions)
    prob = model.problem
    with model:
        # if there's an objective set, make sure there's at least one solution
        if model.objective != Objective(Zero):
            model.slim_optimize(
                error_value=None,
                message="There is no optimal solution for the chosen objective!"
            )
            # Add the previous objective as a variable to the model then set it
            # to zero. This also uses the fraction to create the lower/upper
            # bound for the old objective.
            if model.solver.objective.direction == "max":
                fva_old_objective = prob.Variable(
                    "fva_old_objective",
                    lb=fraction_of_optimum * model.solver.objective.value,
                )
            else:
                fva_old_objective = prob.Variable(
                    "fva_old_objective",
                    ub=fraction_of_optimum * model.solver.objective.value,
                )
            fva_old_obj_constraint = prob.Constraint(
                model.solver.objective.expression - fva_old_objective,
                lb=0,
                ub=0,
                name="fva_old_objective_constraint",
            )
            model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])
        # if pfba_factor was passed, do pFBA to get min total flux, multiply by
        # pfba_factor, and add constraint requiring that total flux not be above
        # that number
        if pfba_factor is not None:
            if pfba_factor < 1.0:
                warn(
                    "The 'pfba_factor' should be larger or equal to 1.",
                    UserWarning,
                )
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum,
                    lb=0,
                    ub=0,
                    name="flux_sum_constraint",
                )
            model.add_cons_vars([flux_sum, flux_sum_constraint])
        # make sure model has no objective function
        model.objective = Zero
        # prepare output dataframe so that it can be updated with each max or
        # min flux as it's computed
        fva_result = pd.DataFrame({
            'minimum' : np.zeros(num_reactions, dtype = float),
            'maximum' : np.zeros(num_reactions, dtype = float)
        }, index = reaction_ids)
        # get minimum fluxes and maximum fluxes in two separate batches so that
        # you can share the same model object between all processes rather than
        # making a new copy of the model for each process so they can all
        # independently set the objective direction
        for direction in ("minimum", "maximum"):
            # unlike the original implementation, use a ProcessPool even if
            # there's only one thread to handle timeouts
            pebble_fva(
                model, reaction_ids, direction, fva_result, processes, verbose,
                update_freq
            )
            # pebble_fva updates the fva_result dataframe directly
    # add a column with the reaction string unless any of the given reaction
    # "IDs" were actually strings of multiple IDs separated by +s
    if not any('+' in r for r in reaction_ids):
        fva_result['reaction'] = fva_result.index.map(
            lambda r: model.reactions.get_by_id(r).build_reaction_string(True)
        )
    return(fva_result)

def pebble_fva(model, reactions, direc, out_df, threads, verbose, update_freq):
    '''
    Set up pebble.ProcessPool to find max or min fluxes for all given reactions.

    When predicting fluxes from a model a bunch of times in a row, sometimes
    the process just hangs indefinitely. The original Cobrapy FVA implementation
    used multiprocessing to minimize or maximize the fluxes through multiple
    reactions in parallel, and multiprocessing does not let you kill a single
    hanging process in a pool of multiple processes that may have several
    totally fine processes, so effectively the whole pool winds up hanging.
    Pebble lets you set a timeout for each individual process, and if any
    process is still running after that amount of time, it's killed and then
    replaced with a new one from the original list. Making a list of all the
    killed processes and then redoing them in a second ProcessPool almost always
    gets everything on the second try, so whatever the issue is appears to be
    entirely independent of anything about the input
    '''
    model.solver.objective.direction = direc[:3]
    # initialize the ProcessPool
    pool = ProcessPool(
        max_workers = threads,
        initializer = _pool_init,
        initargs = (model,)
    )
    to_do = reactions.copy()
    attempts = 0
    max_attempts = 3
    while to_do:
        attempts += 1
        batch_start = time.time()
        if verbose > 1:
            msg = f'Attempting to find {direc} fluxes for {len(to_do)} '
            msg += 'reactions...'
            print(msg)
        # start with a timeout of 5 minutes then up it to 10 and then 30 on
        # subsequent attempts to make sure we don't waste time on hanging
        # processes but also don't cut off the few reactions that genuinely
        # need a while to handle
        if attempts == 1:
            timeout = 300
        elif attempts == 2:
            timeout = 600
        elif attempts == 3:
            timeout = 1800
        future = pool.map(_fva_step, to_do, timeout = timeout)
        iterator = future.result()
        # results will be returned in the order they were submitted in, so keep
        # an iteration counter to know which reactions to redo
        i = 0
        to_redo = list()
        done = 0
        failed = 0
        while True:
            try:
                # no errors; update output dataframe
                done += 1
                out_df.at[to_do[i], direc] = next(iterator)
            except (StopIteration, IndexError):
                # exit the loop if every reaction was at least attempted
                break
            except TimeoutError as error:
                # keep track of which reactions timed out
                rxn = to_do[i]
                to_redo.append(rxn)
                if verbose > 1:
                    msg = f'Took over {timeout/60} minutes to find {direc} flux'
                    msg += f' for {rxn}; will try again unless this is the last'
                    msg += ' attempt'
                    print(msg)
            except ProcessExpired as error:
                # keep track of which reactions encountered other errors
                # without stopping the whole thing
                failed += 1
                if verbose > 1:
                    rxn = to_do[i]
                    msg = f'Couldn\'t get {direc} flux for {rxn}.'
                    print(msg)
                    print(error)
            except Exception as error:
                # keep track of which reactions encountered other errors
                # without stopping the whole thing
                failed += 1
                if verbose > 1:
                    rxn = to_do[i]
                    msg = f'Couldn\'t get {direc} flux for {rxn}.'
                    print(msg)
                    print(error)
            finally:
                # increment iterator regardless of success or failure
                i += 1
                if (i % update_freq == 0) and (verbose > 1):
                    msg = f'On reaction {i} of {len(to_do)} after '
                    msg += f'{round((time.time() - batch_start)/60)} minute(s)'
                    print(msg)
        # all reactions in to_do have now been attempted, but some might need
        # to be redone. If none need redoing, to_do will be empty and the while
        # loop will exit
        if verbose > 1:
            done = len(to_do) - len(to_redo)
            batch_time = round((time.time() - batch_start)/60)
            msg = f'Took {batch_time} minute(s) to attempt to get {direc} '
            msg += f'fluxes for {len(to_do)} reactions. Succeeded for {done}, '
            msg += f'need to retry {len(to_redo)}, and failed for {failed}.'
            print(msg)
        to_do = to_redo.copy()
        if (len(to_do) > 0) and (attempts >= max_attempts):
            if verbose > 0:
                msg = f'After {attempts} tries, did not manage to get {direc} '
                msg += f'for {len(to_redo)} reactions in under {timeout/60} '
                msg += f'minutes; setting {direc} flux for those reactions to NA'
                print(msg)
            for rxn_id in to_redo:
                out_df.at[rxn_id, direc] = float('nan')
            break
    # close the ProcessPool
    pool.close()
    pool.join()
    # don't have to return out_df because we modified it directly

def _pool_init(m):
    '''
    Initialize the Cobrapy Model object as a global variable to be shared by all
    processes in the ProcessPool so we don't waste a bunch of time unnecessarily
    copying it around to each one
    '''
    global model
    model = m

def _fva_step(reaction_id):
    '''
    Compute the maximum or minimum possible flux for the reaction with the
    given ID. Expects the Cobrapy Model object that contains the referenced
    reaction to be a global variable called "model".
    Can also accept strings of reaction IDs separated by +s and/or -s that
    specify a linear combination of reactions to find the maximum or minimum
    net flux through
    '''
    # interpret +s in reaction IDs to be delimiting a list of reaction IDs that
    # we want to find the min or max net flux through
    coefs = dict()
    for r_id in reaction_id.split('+'):
        # if any reaction ID starts with a -, assume that's not actually part of
        # the reaction ID but an indication that the flux through that reaction
        # should be subtracted from the other reactions' fluxes (not added)
        if not r_id.startswith('-'):
            rxn = model.reactions.get_by_id(r_id)
            coefs[rxn.forward_variable] = 1
            coefs[rxn.reverse_variable] = -1
        else:
            rxn = model.reactions.get_by_id(r_id.lstrip('-'))
            coefs[rxn.forward_variable] = -1
            coefs[rxn.reverse_variable] = 1
    model.solver.objective.set_linear_coefficients(coefs)
    model.slim_optimize()
    # make sure solver status was optimal
    if model.solver.status == OPTIMAL:
        value = model.solver.objective.value
    elif model.solver.status == UNBOUNDED:
        # this means the objective value was infinite
        if model.solver.objective.direction == 'max':
            value = float('inf')
        elif model.solver.objective.direction == 'min':
            value = -float('inf')
    else:
        # if it wasn't optimal or unbounded, use NaN
        value = float('nan')
    # revert all objective coefficients we changed at the beginning
    zeros = dict()
    for r_id in reaction_id.split('+'):
        rxn = model.reactions.get_by_id(r_id.lstrip('-'))
        zeros[rxn.forward_variable] = 0
        zeros[rxn.reverse_variable] = 0
    model.solver.objective.set_linear_coefficients(zeros)
    return(value)
