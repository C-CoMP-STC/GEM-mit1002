# modelmaker_dilution.py
'''
Adds a "dilution reaction" for each metabolite in the model (or a given list of
metabolites in the given model) and constrains the flux through that reaction to
be 1/dil_factor * the sum of the absolute values of the fluxes through all
reactions that involve that metabolite weighted by that metabolite's
stoichiometric coefficients in each reaction.

Requires the model to be capable of net production of all metabolites that
participate in at least one reaction with nonzero predicted flux; this prevents
"perfect" recycling of cofactors by requiring at least some net production of
them. This can improve predictions of knockout phenotypes by requiring fluxes
through biosynthetic pathways that might otherwise not be required to have flux.

Mitigates but does not completely solve the problem of unbounded loop fluxes.

This concept was introduced in https://doi.org/10.1186/gb-2010-11-4-r43 and this
implementation of the idea is based on the implementation presented in
https://doi.org/10.1371/journal.pcbi.1003126
'''

import cobra
import re
import pandas as pd
import itertools as it
from optlang.symbolics import Zero

def add_dilution_reactions(
    given_model, mets_to_dilute = None, add_leakage = True, leak_flux = 1,
    verbose = 1
):
    '''
    Given a Cobrapy Model object and, optionally, a list of IDs of Metabolite
    objects in that Model, add a "dilution" reaction for each metabolite in
    mets_to_dilute that irreversibly consumes that metabolite and produces
    nothing to represent "dilution" of that metabolite out of/away from the
    cell, or loss to a side reaction (e.g. oxidative damage).

    If no mets_to_dilute are given, will create one for all metabolites in the
    model except for ones that have "tRNA" or "cytochrome" (case-insensitive)
    in their IDs and/or names.
        Those metabolites are frequently present in models but only capable
        of participating in short loops of reactions (e.g. redox loops for
        cytochromes, charging with an amino acid then losing the amino acid in
        a "translation" reaction for tRNAs) and lack biosynthesis pathways. If
        add_dilution_constraints() is subsequently called on this model, those
        dilution constraints will prevent all metabolites that have dilution
        reactions from being infinitely recycled by such loops without any
        contributing biosynthetic flux. Since tRNA and cytochrome metabolites
        tend to lack biosynthetic pathways in models, cytocrhomes tend to be
        involved in particularly important reactions (e.g. electron transport
        chain), and flux through tRNA-dependent reactions tends to be required
        for biomass flux, it's usually best to just skip them.
    add_leakage governs whether or not add_leakage_reactions is also called; see
    that docstring for details.
    Modifies a copy of the given model and returns the copy
    '''
    # start by making sure given_model doesn't already have dilution reactions
    probably_dilution_rxns = [
        r for r in given_model.reactions
        if r.id.endswith('_dilution') and (len(r.metabolites) == 1)
    ]
    if any(probably_dilution_rxns):
        if verbose > 0:
            msg = 'The given_model passed to add_dilution_reactions appears to '
            msg += 'already have at least one dilution reaction; returning the '
            msg += 'model as-is'
            print(msg)
        return(given_model)
    # otherwise, make a copy of given_model to add dilution reactions to
    model = given_model.copy()
    # if no list of metabolites to add dilution constraints for was given, use
    # the list of all metabolites in the model except for ones with "tRNA" in
    # their names, since those usually lack a biosynthetic pathway (since it's
    # extremely complicated and frequently different for each tRNA) and can thus
    # only participate in a loop between being charged with an amino acid and
    # being returned to their original state by a "translation" reaction, which
    # tends to be required for biomass flux
    if mets_to_dilute is None:
        to_skip = ['trna', 'cytochrome']
        mets_to_dilute = [
            m.id for m in model.metabolites
            # be case-insensitive
            if all(s not in m.name.lower() for s in to_skip)
            # check metabolite names and IDs, just to be thorough
            and all(s not in m.id.lower() for s in to_skip)
        ]
        if verbose > 0:
            msg = 'Since no mets_to_dilute were passed to '
            msg += 'add_dilution_reactions, adding dilution reactions for all '
            msg += f'{len(mets_to_dilute)} metabolites that don\'t have "tRNA" '
            msg += 'or "cytochrome" in their names or IDs. See docstring for '
            msg += 'explanation.'
            print(msg)
    elif not all(isinstance(m, str) for m in mets_to_dilute):
        # if this was a list of Metabolite objects from given_model, replace it
        # with a list of their IDs, since those objects are technically not
        # the same as the metabolite objects in model, and that'll cause drama
        # later
        try:
            mets_to_dilute = [m.id for m in mets_to_dilute]
        except:
            msg = 'Something was wrong with the contents of mets_to_dilute in '
            msg += 'add_dilution_reactions(); make sure it is either a list of '
            msg += 'cobra.Metabolite objects from given_model or a list of '
            msg += 'strings corresponding to the IDs of cobra.Metabolite '
            msg += 'objects that are present in given_model.'
            raise ValueError(msg)
    # now make an irreversible reaction that consumes each metabolite in
    # mets_to_dilute whose ID is <metabolite ID>_dilution so they're easy to
    # identify later
    dil_rxns = list()
    for met_id in mets_to_dilute:
        met_obj = model.metabolites.get_by_id(met_id)
        dil_rxn = cobra.Reaction(f'{met_id}_dilution')
        dil_rxn.name = f'{met_obj.name} Dilution'
        dil_rxn.lower_bound = 0
        dil_rxn.upper_bound = float('Inf')
        dil_rxn.add_metabolites({met_obj : -1.0})
        dil_rxns.append(dil_rxn)
    # adding a list of many reactions is much faster than separately adding
    # many individual reactions due to the way Cobrapy Models work
    model.add_reactions(dil_rxns)
    # add "leakage" reactions if requested
    if add_leakage:
        model = add_leakage_reactions(model, leak_flux, verbose)
    return(model)

def add_leakage_reactions(given_model, bound, verbose):
    '''
    Finds all pairs of Metabolite objects in given_model that represent the
    same real-world compound in different subcellular compartments and creates
    a new reversible reaction that reversibly interconverts them with a flux of
    up to +/- bound (which should be a relatively small non-negative non-zero 
    number). This represents "leakage" of this metabolite across the boundary
    separating the two compartments.

    This is important if the only way for a particular pair of metabolites to
    to move between a particular pair of compartments is via antiport, i.e.
    A[out] + B [in] <-> A[in] + B[out]
    and the only way to produce B[in] is to derive it from A[out], which might
    be the case for an NTP and its corresponding NDP.

    Dilution constraints would require a small amount of A[in] to be "wasted"
    via its dilution reaction, thus ensuring that there is always less B[in]
    than the amount of A[out] that initially entered, so there is never enough
    B[in] to exchange for more A[out], preventing the antiport reaction from 
    ever having fluxes greater than 0.

    Leakage reactions resolve this conundrum this by allowing an arbitrary
    small flux of each metabolite to move between compartments independently
    of any possible antiport scheme to make up for these dilution fluxes. They
    are also hypothetically plausible representations of reality, since no
    biological membrane is always completely impemetrable to any particular
    metabolite, and we're usually trying to model steady states of metabolic
    networks.

    A good choice for the value of bound is 1/the dilution_factor you pass to
    add_dilution_constraints()
    '''
    # first, make sure given_model doesn't already have leakage reactions
    if any(
        ('--' in r.id) and r.id.endswith('_leakage')
        for r in given_model.reactions
    ):
        if verbose > 0:
            msg = 'The given_model passed to add_leakage_reactions appears to '
            msg += 'already contain leakage reactions; returning the model '
            msg += 'as-is.'
            print(msg)
        return(given_model)
    # otherwise, as usual, modify a copy of the given model
    model = given_model.copy()
    leakage_rxns = list()
    # in a perfect world, metabolites in different compartments that represent
    # the same real-world metabolite would contain some reference to each other
    # in their IDs, but we do not live in a perfect world, so hope that these
    # metabolites can be linked with their names. start by making a list of all
    # metabolites that have the same name as at least one other metabolite
    all_names = list()
    for m in model.metabolites:
        # if metabolite names appear to end with their compartments, drop the
        # compartment names so we still find matches
        regexp = re.compile(r' ?[\(\[\{]?' + m.compartment + r'[\)\]\}]?$')
        all_names.append(regexp.sub('', m.name))
    name_counts = pd.Series(all_names).value_counts()
    shared_names = name_counts[name_counts > 1].index.unique().tolist()
    # now loop over this list of names and get the corresponding metabolites
    for met_name in shared_names:
        met_objs = [m for m in model.metabolites if met_name == m.name]
        # loop over each possible pair of these metabolites
        for (m1, m2) in it.combinations(met_objs, 2):
            # a Cobrapy Metabolite object has a "reactions" attribute that is a
            # frozenset of all the Cobrapy Reaction objects that involve that
            # Metabolite object
            shared_rxns = m1.reactions.intersection(m2.reactions)
            # see if any reaction involves both of these metabolites
            if shared_rxns:
                # create a reversible reaction that interconverts these
                leak_rxn = cobra.Reaction(
                    f'{m1.id}--{m2.id}_leakage', name = f'{met_name} Leakage'
                )
                leak_rxn.add_metabolites({m1 : -1, m2 : 1})
                # only let it have a flux of up to 1/dil_factor
                leak_rxn.lower_bound = -bound
                leak_rxn.upper_bound = bound
                leakage_rxns.append(leak_rxn)
    if verbose > 0:
        msg = f'Added a "leakage" reaction for all {len(leakage_rxns)} pairs '
        msg += 'of metabolites in different compartments that represent the '
        msg += 'same real-world metabolite and set the bounds on each new '
        msg += f'reaction to +/- {bound}.'
        print(msg)
    # add all leakage reactions in one go cuz that's much faster than adding
    # them one at a time
    model.add_reactions(leakage_rxns)
    return(model)

def add_dilution_constraints(
    given_model, dil_factor = 0, debug = False, debug_rxn = '', verbose = 1
):
    '''
    Given a Cobrapy Model object that contains at least one dilution reaction
    created by add_dilution_reactions and some number larger than 1, create a
    constraint for each dilution reaction that requires the flux through that
    dilution reaction to be equal to the sum of the fluxes through all other
    reactions that involve that metabolite divided by dil_factor.
    Modifies a copy of the given model and returns the copy
    '''
    # first make sure the model has at least one dilution reaction to create a
    # constraint for
    diluted_mets = [
        list(r.metabolites)[0].id for r in given_model.reactions
        if r.id.endswith('_dilution') and (len(r.metabolites) == 1)
    ]
    if len(diluted_mets) == 0:
        msg = 'The model passed to add_dilution_constraints does not appear to '
        msg += 'have any dilution reactions in it already; try calling '
        msg += 'add_dilution_reactions() on this model first.'
        raise ValueError(msg)
    # then make sure dil_factor is > 0, cuz there's no point in making dilution
    # constraints with dilution factors of 0 (they'd divide by zero)
    if dil_factor <= 0:
        if verbose > 0:
            msg = 'The dil_factor passed to add_dilution_constraints was 0; '
            msg += 'returning given_model without dilution constraints.'
            print(msg)
        return(given_model)
    # otherwise, modify a copy of the given model
    model = given_model.copy()
    if debug:
        # if debug is True, add each dilution reaction one at a time,
        # check to see if the given reporter reaction can still
        # sustain flux, and stop if/when a particular dilution
        # constraint blocks flux through the given reporter reaction
        model = debug_dilution(
            model, debug_rxn, diluted_mets, dil_factor
        )
    else:
        # otherwise, just make all the dilution constraints in a list comp
        dil_consts = [
            _make_dilution_constraint(model, met_id, dil_factor)
            for met_id in diluted_mets
        ]
        model.add_cons_vars(dil_consts)
    return(model)

def _make_dilution_constraint(model, met_id, dil_factor):
    # start out with an optlang/SymPy expression that's just zero, and add in
    # the optlang/SymPy variables associated with the reactions that this
    # metabolite participates in
    expression = Zero
    # passing around Metabolite objects directly has frequently given us weird
    # errors, so we're passing around the metabolite IDs instead
    met_obj = model.metabolites.get_by_id(met_id)
    for r in met_obj.reactions:
        if 'dilution' not in r.id:
            # add both the forward and reverse variable so this is always
            # positive, regardless of which direction the reaction is going in
            expression += r.forward_variable + r.reverse_variable
        else:
            # if this is the metabolite's dilution reaction, subtract its
            # flux times the dilution factor (since it's irreversible, reverse
            # variable should always be 0, but subtracted just in case)
            expression -= dil_factor * (r.forward_variable - r.reverse_variable)
    # set the upper and lower bounds on this constraint to 0 so that the
    # dilution flux (scaled by the dilution factor) must equal the sum of fluxes
    # through all other reactions involving this metabolite
    dilution_constraint = model.problem.Constraint(
        expression, lb = 0, ub = 0, name = f'{met_id}_dilution_constraint'
    )
    return(dilution_constraint)

def debug_dilution(model, ref_rxn, mets_to_dilute, dil_factor):
    '''
    Use a pebble ProcessPool with only one process to add each dilution
    constraint, one at a time, and see when the given reaction can no longer
    sustain flux. The ProcessPool is for catching when solving the model hangs
    after imposing a particular set of dilution constraints, which seems to
    happen at some point every single time I run this on models with more than
    a few hundred metabolites. No idea why.
    '''
    import random
    import pebble as pb
    from concurrent.futures import TimeoutError
    # set ref_rxn as the model's objective
    rxn_obj = model.reactions.get_by_id(ref_rxn)
    model.solver.objective.set_linear_coefficients({
        rxn_obj.forward_variable : 1, rxn_obj.reverse_variable : -1
    })
    msg = 'Adding dilution constraints one at a time and doing FBA in between '
    msg += f'adding each constraint to see when {ref_rxn}: '
    msg += f'{rxn_obj.build_reaction_string(True)} can no longer sustain flux.'
    print(msg)
    # randomize the list of metabolites to dilute so we can tell if a particular
    # constraint always reduces flux through ref_rxn to 0 or if several can by
    # reruning this multiple times
    random.shuffle(mets_to_dilute)
    # make a list of lists where the first list has the first metabolite ID,
    # the second has the first two, the third has the first three, etc.
    args = [
        (model, mets_to_dilute[:i], dil_factor)
        for i in range(len(mets_to_dilute))
    ]
    # make a process pool with the model and reaction ID as global objects
    # since they'll stay the same in each iteration
    pool = pb.ProcessPool(max_workers = 1)
    # if it takes more than 30 seconds to find the maximum flux after imposing
    # a particular set of dilution constraints, assume it's infeasible
    future = pool.map(_test_dil_const, args, timeout = 30)
    iterator = future.result()
    failed = False
    i = 0
    while True:
        try:
            # if there were no errors, see if the maximum possible flux is > 0
            max_flux = next(iterator)
            if abs(max_flux) < 10 ** -8:
                met_id = mets_to_dilute[i]
                met_obj = model.metabolites.get_by_id(met_id)
                met_str = f'{met_obj.name} [{met_obj.compartment}]'
                rxn_obj = model.reactions.get_by_id(ref_rxn)
                rxn_str = rxn_obj.build_reaction_string(True)
                msg = f'Flux through {ref_rxn}: {rxn_str} dropped to 0 after '
                msg += f'imposing dilution constraint for {met_str} ({met_id}) '
                msg += f'after imposing {i} of {len(mets_to_dilute)} dilution '
                msg += 'constraints.'
                print(msg)
                failed = True
                break
        except StopIteration:
            # should only happen if we've imposed all dilution constraints
            break
        except TimeoutError:
            met_id = mets_to_dilute[i]
            met_obj = model.metabolites.get_by_id(met_id)
            met_str = f'{met_obj.name} [{met_obj.compartment}]'
            rxn_obj = model.reactions.get_by_id(ref_rxn)
            rxn_str = rxn_obj.build_reaction_string(True)
            msg = f'Timed out when finding maximum flux through {ref_rxn}: '
            msg += f'{rxn_str} after imposing dilution constraint for '
            msg += f'{met_str} ({met_id}) after imposing {i} of '
            msg += f'{len(mets_to_dilute)} dilution constraints.'
            print(msg)
        except pb.ProcessExpired as error:
            met_id = mets_to_dilute[i]
            met_obj = model.metabolites.get_by_id(met_id)
            met_str = f'{met_obj.name} [{met_obj.compartment}]'
            rxn_obj = model.reactions.get_by_id(ref_rxn)
            rxn_str = rxn_obj.build_reaction_string(True)
            msg = f'Failed to find maximum flux through {ref_rxn}: {rxn_str} '
            msg += f'after imposing dilution constraint for {met_str} '
            msg += f'({met_id}) after imposing {i} of {len(mets_to_dilute)} '
            msg += 'dilution constraints.'
            print(msg)
            failed = True
            break
        except Exception as error:
            met_id = mets_to_dilute[i]
            met_obj = model.metabolites.get_by_id(met_id)
            met_str = f'{met_obj.name} [{met_obj.compartment}]'
            rxn_obj = model.reactions.get_by_id(ref_rxn)
            rxn_str = rxn_obj.build_reaction_string(True)
            msg = f'Failed to find maximum flux through {ref_rxn}: {rxn_str} '
            msg += f'after imposing dilution constraint for {met_str} '
            msg += f'({met_id}) after imposing {i} of {len(mets_to_dilute)} '
            msg += 'dilution constraints.'
            print(msg)
            print(error)
            print(error.traceback)
            failed = True
            break
        finally:
            i += 1
            if (i % 1) == 0:
                print(f'On constraint {i} of {len(mets_to_dilute)}')
    if failed:
        raise Exception
    else:
        rxn_obj = model.reactions.get_by_id(ref_rxn)
        rxn_str = rxn_obj.build_reaction_string(True)
        msg = f'Managed to successfully impose all {len(dil_consts)} '
        msg += 'dilution constraints while maintaining non-zero flux through '
        msg += f'{rxn}: {rxn_str}'
        print(msg)
        model.add_cons_vars(dil_consts)
        return(model)

def _test_dil_const(args):
    (model, met_ids, dil_factor) = args
    model.add_cons_vars([
        make_dilution_constraint(model, m, dil_factor) for m in met_ids
    ])
    max_flux = model.slim_optimize()
    return(max_flux)
