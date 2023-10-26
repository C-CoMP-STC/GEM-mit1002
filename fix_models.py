# fix_models.py
'''
Functions used by the different fix_<name of generic GSMM>.py scripts
'''

def fix_diphosphate_rxns(model, ppi_ids, pi_ids, verbose):
    '''
    ppi_ids is the list of metabolite IDs that correspond to diphosphate ions
    pi_ids is the list of metabolite IDs that correspond to phosphate ions
    Nominally, hydrolysis of a molecule with a diphosphate group is almost
    perfectly reversible at physiological pHs, but real cells express lots of
    highly active diphosphatases specifically to drive this reaction in the
    direction of hydrolysis. This is crucial in real life because DNA and RNA
    polymerases catalyze the reaction
    DNA/RNA + (d)NTP + H2O <-> DNA/RNA-NMP + PPi
    and crucial for steady-state modeling, since, unless you're explicitly
    including thermodynamic constraints, things like FBA will abuse the
    reversibility of these reactions to generate ATP from dozens of utterly
    implausible sources
    '''
    msgs = list()
    for r in model.reactions:
        # skip reactions that involve diphosphate and phosphate cuz they're
        # probably diphosphatase reactions
        pi_in_rxn = any(m.id in pi_ids for m in r.metabolites)
        # also skip exchange and irreversible reactions
        if (r not in model.boundary) and not pi_in_rxn and r.reversibility:
            # figure out if diphosphate is a product or reactant (or both)
            ppi_prod = any(m.id in ppi_ids for m in r.products)
            ppi_reac = any(m.id in ppi_ids for m in r.reactants)
            if ppi_prod and not ppi_reac:
                r.lower_bound = 0
                msg = f'Irreversibilized diphosphate reaction {r.id}: '
                msg += r.build_reaction_string(True)
                msgs.append(msg)
            elif ppi_reac and not ppi_prod:
                # we could just set the upper bound to 0 and call it a day, but
                # then these reactions would always have negative (or 0) fluxes,
                # and that's a little annoying, so swap the products and
                # reactants and then set the lower bound to 0
                new_met_dict = {
                    # negate all coefficients to swap prods & reacts
                    met : -1 * stoich for (met, stoich) in r.metabolites.items()
                }
                # you can only "add" metabolites to reactions, so do this twice
                # (once will just make their coefficients 0)
                r.add_metabolites(new_met_dict)
                r.add_metabolites(new_met_dict)
                # now we can set the lower bound to 0
                r.lower_bound = 0
                msg = 'Flipped & irreversibilized diphosphate reaction '
                msg += f'{r.id}: {r.build_reaction_string(True)}'
                msgs.append(msg)
    if verbose:
        # sort these messages by reaction ID cuz there may be quite a few
        msgs.sort()
        for msg in msgs:
            print(msg)
    # no particular need to return the model since it was modified in-place

def find_dead_ends(model, verbose):
    '''
    Find all metabolites that have 0 or 1 reactions associated with them, as
    well as any reactions involving those metabolites
    '''
    rxns_to_remove = list()

    # loop over metabolites and assess their dead-ended-ness and the dead-ended
    # ness of any reactions that involve them
    mets_to_remove = list()
    for met in model.metabolites:
        _find_dead_ends_inner(met, mets_to_remove, rxns_to_remove, verbose)
    model.remove_reactions(rxns_to_remove)
    model.remove_metabolites(mets_to_remove)
    if verbose > 0:
        msg = f'Removed {len(rxns_to_remove)} reactions and '
        msg += f'{len(mets_to_remove)} metabolites in total for being dead-ends'
        print(msg)

def _find_dead_ends_inner(met_to_check, isolated_mets, dead_ends, verbose):
    '''
    Check if a given metabolite is only associated with a single reaction, then
    check all the other metabolites in that reaction but avoid checking the
    same metabolite twice so it doesn't recurse infinitely
    '''
    # first make sure this isn't already on our list of isolated metabolites
    if met_to_check not in isolated_mets:
        # if any reactions this metabolite participates in are already on our
        # list of dead-ends, don't count them when deciding if this metabolite
        # is also a dead-end/isolated
        rxns_to_check = [
            r for r in met_to_check.reactions if r not in dead_ends
        ]
        # now see if the metabolite still has at least 2 non-dead-end reactions
        if len(rxns_to_check) < 2:
            # this metabolite is only associated with 0 or 1 reactions, so it's
            # a dead-end/isolated metabolite
            isolated_mets.append(met_to_check)
            if verbose:
                msg = f'Removed isolated/dead-end metabolite {met_to_check.id}:'
                msg += f' {met_to_check.name}'
                print(msg)
            # if there's only one reaction associated with this metabolite, it's
            # a dead-end
            if len(rxns_to_check) == 1:
                dead_ends.append(rxns_to_check[0])
                if verbose:
                    msg = f'Removed dead-end reaction {rxns_to_check[0].id}: '
                    msg += rxns_to_check[0].build_reaction_string(True)
                    print(msg)
                # check the other metabolites in this reaction
                for met in rxns_to_check[0].metabolites:
                    # might as well skip metabolites we already know are
                    # isolated here
                    if met not in isolated_mets:
                        _find_dead_ends_inner(
                            met, isolated_mets, dead_ends, verbose
                        )
        else:
            # this metabolite participates in at least 2 reactions that we
            # haven't already flagged as being dead-ends, but see if they all
            # irreversibly produce or irreversibly consume this metabolite; if
            # so, then they're all dead end reactions
            all_irrev = all(not r.reversibility for r in rxns_to_check)
            all_prod = all(
                met_to_check in r.products for r in rxns_to_check
            )
            all_cons = all(
                met_to_check in r.reactants for r in rxns_to_check
            )
            if all_irrev and (all_prod or all_cons):
                # this metabolite can only ever be produced and never
                # consumed or vice versa, so it and all reactions it
                # participates in are dead-ends
                isolated_mets.append(met_to_check)
                if verbose:
                    msg = f'Removed dead-end metabolite {met_to_check.id}: '
                    msg += f'{met_to_check.name}'
                    print(msg)
                for r in rxns_to_check:
                    dead_ends.append(r)
                    if verbose:
                        msg = f'Removed dead-end reaction {r.id}: '
                        msg  += r.build_reaction_string(True)
                        print(msg)
                # check all the other metabolites that participate in these
                # reactions to see if removing these reactions would also
                # isolate them
                mets_to_check = set([
                    m for r in rxns_to_check
                    for m in r.metabolites
                    if m not in isolated_mets
                ])
                for m in mets_to_check:
                    _find_dead_ends_inner(
                        m, isolated_mets, dead_ends, verbose
                    )

def find_secretly_irrev_rxns(model):
    '''
    Some reactions have both their upper and lower bounds set to non-zero values
    but have at least one metabolite that can only be consumed or only be
    produced by irreversible reactions. This function identifies all such
    reactions and sets their lower bounds to 0 if they're only capable of going
    forwards or swaps their products and reactants and then sets their lower
    bounds to zero if they were originally only capable of going backwards
    '''
    msgs = list()
    i = 0
    for rev_rxn in model.reactions:
        # skip exchange reactions just because
        if rev_rxn.reversibility and not rev_rxn.boundary:
            # once we've found one metabolite that can only be consumed or only
            # be produced, we don't need to check the rest
            changed = False
            for met in rev_rxn.metabolites:
                msg = _find_secretly_irrev_rxns_inner(met, rev_rxn)
                msgs.append(msg)
                # msg will be an empty string unless we made this reaction
                # irreversible
                if msg != '':
                    changed = True
                if changed:
                    i += 1
                    break
    # sort messages before printing them so it's less chaotic
    msgs.sort()
    for msg in msgs:
        if msg != '':
            print(msg)
    msg = f'In total, made {i} reversible reactions irreversible because they '
    msg += 'had at least one product that could only be produced by other '
    msg += 'reactions or one reactant that could only be consumed by other '
    msg += 'reactions.'
    print(msg)
    # everything was modified in-place, so no need to return anything

def _find_secretly_irrev_rxns_inner(met, rev_rxn):
    # make a list of all other reactions this metabolite participates in
    other_rxns = [r for r in met.reactions if r != rev_rxn]
    # see if all of those reactions are irreversible and either all consume or
    # all produce this metabolite
    all_irrev = all(not r.reversibility for r in other_rxns)
    all_prod = all(met in r.products for r in other_rxns)
    all_cons = all(met in r.reactants for r in other_rxns)
    if all_irrev and all_cons and (met in rev_rxn.products):
        # A <-> B ->
        rxn_str = rev_rxn.build_reaction_string(True)
        rev_rxn.lower_bound = 0
        msg = f'Set lower bound to 0 for {rev_rxn.id}: {rxn_str} because all '
        msg += f'other reactions involving {met.name} ({met.id}) could only '
        msg += 'consume it and never produce it'
    elif all_irrev and all_prod and (met in rev_rxn.reactants):
        # -> A <-> B
        rxn_str = rev_rxn.build_reaction_string(True)
        rev_rxn.lower_bound = 0
        msg = f'Set lower bound to 0 for {rev_rxn.id}: {rxn_str} because all '
        msg += f'other reactions involving {met.name} ({met.id}) could only '
        msg += 'produce it and never consume it'
    elif all_irrev and all_cons and (met in rev_rxn.reactants):
        # <- A <-> B
        rxn_str = rev_rxn.build_reaction_string(True)
        # swap products and reactants before setting lower bound to 0
        new_met_dict = {
            # multiplying current stoichiometric coefficients by -2 and adding
            # to existing ones cuz we apparently can only add or subtract
            m : -2 * stoich for (m, stoich) in rev_rxn.metabolites.items()
        }
        rev_rxn.add_metabolites(new_met_dict)
        rev_rxn.lower_bound = 0
        msg = 'Switched products and reactants and set lower bound to 0 for '
        msg += f'{rev_rxn.id}: {rxn_str} because all other reactions involving '
        msg += f'{met.name} ({met.id}) could only consume it and never produce '
        msg += 'it'
    elif all_irrev and all_prod and (met in rev_rxn.products):
        # A <-> B <-
        rxn_str = rev_rxn.build_reaction_string(True)
        new_met_dict = {
            m : -2 * stoich for (m, stoich) in rev_rxn.metabolites.items()
        }
        rev_rxn.add_metabolites(new_met_dict)
        rev_rxn.lower_bound = 0
        msg = 'Switched products and reactants and set lower bound to 0 for '
        msg += f'{rev_rxn.id}: {rxn_str} because all other reactions involving '
        msg += f'{met.name} ({met.id}) could only produce it and never consume '
        msg += 'it'
    else:
        msg = ''
    # reaction was modified in-place, but return the message so we can sort them
    # by type before printing so they're a bit more organized
    return(msg)

def drop_stoich_gpr_dupes(model, to_drop, verbose):
    '''
    Finds pairs of reactions that involve exactly the same metabolites but with
    different stoichiometry, reversibility, and/or GPRs
    Merges each pair into a single reaction in various ways (see output to see
    exactly how and why it changes what it does)
    '''
    # loop over all pairs of reactions
    better_rxns = list()
    for (r1, r2) in it.combinations(model.reactions, 2):
        # check if they involve exactly the same metabolites
        if set(r1.metabolites.keys()) == set(r2.metabolites.keys()):
            # to check if all the stoichiometric coefficients in one reaction
            # are n * the corresponding stoichiometric coefficient in the other
            # for any n, make a matrix with one column for each reaction and one
            # row for each metabolite and see what the rank of that matrix is
            stoich_mat = np.array([
                (r1.metabolites[met], r2.metabolites[met])
                # make sure both columns are in the same order
                for met in r1.metabolites.keys()
            ])
            if np.linalg.matrix_rank(stoich_mat) == 1:
                # see if all stoichiometric coefficients are the same magnitude
                # but opposite sign
                sign_diff_only = all(
                    r1.metabolites[met] == -r2.metabolites[met]
                    for met in r1.metabolites.keys()
                )
                # next steps depend on the reversibilities of these reactions
                if (not r1.reversibility) and (not r2.reversibility):
                    # if they're both irreversible, see if the only difference
                    # is that they're oriented in opposite directions
                    sign_diff_only = all(
                        r1.metabolites[met] == -r2.metabolites[met]
                        for met in r1.metabolites.keys()
                    )
                    # only need to do something if there are more differences
                    if not sign_diff_only:
                        # I happen to know all of these cases have identical
                        # GPRs but one set of stoichiometric coefficients with
                        # a greatest common divisor that's larger than 1
                        r1_stoichs = [
                            int(abs(x)) for x in r1.metabolites.values()
                        ]
                        #r1_gcd = gcd(*r1_stoichs)
                        r1_gcd = reduce(gcd, r1_stoichs)
                        r2_stoichs = [
                            int(abs(x)) for x in r2.metabolites.values()
                        ]
                        #r2_gcd = gcd(*r2_stoichs)
                        r2_gcd = reduce(gcd, r2_stoichs)
                        # drop the reaction with the GCD that isn't 1
                        if r1_gcd > 1:
                            to_drop.append(r1)
                            better_rxns.append(r2)
                        elif r2_gcd > 1:
                            to_drop.append(r2)
                            better_rxns.append(r1)
                elif r1.reversibility and r2.reversibility:
                    # if they're both reversible, see if they have the same GPR
                    if r1.gene_reaction_rule == r2.gene_reaction_rule:
                        # I happen to know that if they have the same GPR, they
                        # always have the same coefficients but are oriented in
                        # opposite directions, so arbitrarily pick one to drop
                        r1_num = int(r1.id.lstrip('MAR'))
                        r2_num = int(r2.id.lstrip('MAR'))
                        # drop the one with the higher ID number
                        if r1_num > r2_num:
                            to_drop.append(r1)
                            better_rxns.append(r2)
                        else:
                            to_drop.append(r2)
                            better_rxns.append(r1)
                    else:
                        # drop the reaction associated with fewer genes
                        if len(r1.genes) < len(r2.genes):
                            to_drop.append(r1)
                            better_rxns.append(r2)
                            # divide coefficients by their greatest common
                            # divisor just in case that isn't already 1
                            norm_coefs(r2)
                        else:
                            to_drop.append(r2)
                            better_rxns.append(r1)
                            norm_coefs(r1)
                else:
                    # if one is reversible and one is irreversible, I happen to
                    # know that it's always the case that the irreversible one
                    # is the one that should be deleted
                    if r1.reversibility:
                        to_drop.append(r2)
                        better_rxns.append(r1)
                    elif r2.reversibility:
                        to_drop.append(r1)
                        better_rxns.append(r2)
    # print deets about the reactions flagged as needing to be dropped
    for (bad, good) in zip(to_drop, better_rxns):
        msg = f'Removed {bad.id}\n\t{bad.build_reaction_string(True)}\n'
        msg += f'\t{bad.reaction}\n\tGPR: {bad.gene_name_reaction_rule}\n\t'
        msg += f'GPR: {bad.gene_reaction_rule}\nbecause {good.id} exists\n\t'
        msg += f'{good.build_reaction_string(True)}\n\t{good.reaction}\n\t'
        msg += f'GPR: {good.gene_name_reaction_rule}\n\tGPR: '
        msg += f'{good.gene_reaction_rule}\n'
        if verbose:
            print(msg)

def norm_coefs(rxn):
    '''
    Find greatest common divisor of absolute values of all stoichiometric
    coefficients for this reaction and divide all coefficients by that number
    '''
    rxn_gcd = reduce(gcd, [abs(x) for x in rxn.metabolites.values()])
    #rxn_gcd = gcd(*[int(abs(x)) for x in rxn.metabolites.values()])
    rxn.add_metabolites(
        {met : stoich/rxn_gcd for (met, stoich) in rxn.metabolites.items()},
        combine = False
    )
    # don't have to return anything since reaction was modified in-place

def merge_transport_pairs(model, to_drop, verbose):
    '''
    There are a lot of pairs of transport reactions that are suspiciously
    similar where we really only need to keep one
    '''
    bad_pairs = set()
    for ref_met in list(model.metabolites):
        # find transport reactions involving this metabolite
        trans_rxns = [r for r in ref_met.reactions if any(
            other_met.id[:-1] == ref_met.id[:-1]
            for other_met in r.metabolites
        )]
        # loop over pairs of transport reactions
        for (r1, r2) in it.combinations(trans_rxns, 2):
            # see if the two stoichiometry dicts are either identical or identical
            # if you multiply all coefficients by -1 (i.e. same or opposite
            # directions)
            same_dir = r1.metabolites == r2.metabolites
            opp_dir = r1.metabolites == {
                m : -1 * s for (m,s) in r2.metabolites.items()
            }
            # see if either (solved) GPR is a subset of the other
            r1_gcs = rp.solve(r1.gene_reaction_rule)
            r2_gcs = rp.solve(r2.gene_reaction_rule)
            gpr_1_in_2 = all(g1 in r2_gcs for g1 in r1_gcs)
            gpr_2_in_1 = all(g2 in r1_gcs for g2 in r2_gcs)
            if (same_dir or opp_dir) and (gpr_1_in_2 or gpr_2_in_1):
                # order reactions by ID so we don't add the same pair twice in
                # opposite orders
                if int(r1.id[3:]) < int(r2.id[3:]):
                    bad_pairs.add((r1, r2))
                else:
                    bad_pairs.add((r2, r1))
    # loop over the bad pairs so we don't do this twice
    for (r1, r2) in bad_pairs:
        msg = f'{r1.id}: {r1.build_reaction_string(True)}\nGPR: '
        msg += f'{r1.gene_reaction_rule}\n{r2.id}: '
        msg += f'{r2.build_reaction_string(True)}\nGPR: '
        msg += f'{r2.gene_reaction_rule}\n'
        # check if they're oriented in opposite directions
        opp_dir = r1.metabolites == {
            m : -1 * s for (m,s) in r2.metabolites.items()
        }
        # keep the reaction with the longer GPR
        if len(r1.genes) > len(r2.genes):
            to_drop.append(r2)
            msg += f'Keeping {r1.id} and removing {r2.id}\n'
            # if they're in opposite directions, make the other one reversible
            # if it isn't already
            if opp_dir and not r1.reversibility:
                r1.lower_bound = -1000
                msg = msg.replace(' and r', ', making it reversible, and r')
        elif len(r1.genes) < len(r2.genes):
            to_drop.append(r1)
            msg += f'Keeping {r2.id} and removing {r1.id}\n'
            if opp_dir and not r2.reversibility:
                r2.lower_bound = -1000
                msg = msg.replace(' and r', ', making it reversible, and r')
        else:
            # if the GPRs are the same, see if one is already reversible
            if r1.reversibility and not r2.reversibility:
                to_drop.append(r2)
                msg += f'Keeping {r1.id} and removing {r2.id}\n'
            elif r2.reversibility and not r1.reversibility:
                to_drop.append(r1)
                msg += f'Keeping {r2.id} and removing {r1.id}\n'
            else:
                # if both are reversible or irreversible, drop the one with the
                # higher ID
                if int(r1.id[3:]) < int(r2.id[3:]):
                    to_drop.append(r2)
                    msg += f'Keeping {r1.id} and removing {r2.id}\n'
                    # if both were irreversible and point in opposite
                    # directions, make the one we kept reversible
                    if opp_dir and not r1.reversibility:
                        r1.lower_bound = -1000
                        msg = msg.replace(
                            ' and r', ', making it reversible, and r'
                        )
                else:
                    to_drop.append(r1)
                    msg += f'Keeping {r2.id} and removing {r1.id}\n'
                    if opp_dir and not r2.reversibility:
                        r2.lower_bound = -1000
                        msg = msg.replace(
                            ' and r', ', making it reversible, and r'
                        )
        if verbose:
            print(msg)

def drop_blocked_rxns(model, verbose):
    '''
    Remove all reactions that have both bounds set to 0
    '''
    to_remove = list()
    for r in model.reactions:
        if (r.upper_bound == 0) and (r.lower_bound == 0):
            to_remove.append(r)
            if verbose:
                r_str = r.build_reaction_string(True)
                print(f'Removed blocked reaction {r.id}: {r_str}')
    model.remove_reactions(to_remove)
    if verbose:
        msg = f'Removed {len(to_remove)} reactions that had both their upper '
        msg += 'and lower bounds set to 0'
        print(msg)
    # once again, model was modified in-place, so no particular reason to return
