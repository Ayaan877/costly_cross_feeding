import numpy as np
from find_intermediates import get_candidates
from reverse_scope import giveRevScope
from batch_pruning import *
from satisfiability_check import markSatMetsRxns
from prune_check import isCoreProduced
from autonomy_check import verify_autonomy

def remove_core_production(net, removed_core, rxnMat, prodMat, sumRxnVec,
                           nutrientSet, Currency, Core, intermediates):
    '''
    Prunes a network until it can no longer produce removed_core from base
    nutrients, while keeping all other 7 core metabolites and any specified
    intermediates reachable. Reactions are tested in random order and removed
    greedily whenever the protected set (remaining cores + intermediates) stays
    satisfiable. This creates the "hole" in the network that must be filled by
    a cross-feeding partner.

    Returns the pruned reaction array if removed_core is successfully eliminated,
    or None if no pruning could sever its production without disrupting other
    protected metabolites.
    '''

    remaining_cores = [c for c in Core if c != removed_core]
    protected = list(remaining_cores)
    if intermediates is not None:
        protected.extend(intermediates)

    rxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    rxnVec[net] = 1

    while True:
        remaining = np.nonzero(rxnVec)[0]
        removed_any = False
        for rxn in np.random.permutation(remaining):
            if isCoreProduced(rxn, rxnVec, rxnMat, prodMat,
                              sumRxnVec, nutrientSet, Currency, protected):
                rxnVec[rxn] = 0
                removed_any = True
        if not removed_any:
            break
    
    remaining = np.nonzero(rxnVec)[0]
    satMets, _ = markSatMetsRxns(rxnVec, rxnMat, prodMat,
                                 sumRxnVec, nutrientSet, Currency)

    if satMets[removed_core]:
        return None
    
    return remaining


def make_donor_pathway(intermediate, core_target, rxnMat, prodMat, sumRxnVec, Currency):
    '''
    Constructs a minimal pathway that converts a secreted intermediate into
    a target core metabolite, using only the intermediate itself (plus currency
    metabolites) as the nutrient source. Runs reverse scope expansion from
    core_target with the intermediate as the sole nutrient, then prunes to
    a minimal pathway. An additional check confirms the intermediate is actually
    consumed by at least one reaction in the returned pathway.

    Returns a numpy array of reaction indices for the minimal donor pathway,
    or None if no path exists from the intermediate to the core target.
    '''

    try:
        satMets, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec,
                                        [intermediate], Currency, core_target)
    except ValueError:
        return None
    
    # Check if core is reachable from the intermediate
    if np.sum(satRxns) == 0:
        return None
    
    donor_nutrient = [intermediate]

    donor_pathway = randMinNetwork(satRxns, rxnMat, prodMat, sumRxnVec,
                             core_target, donor_nutrient, Currency)
    
    # Failsafe against pathways that don't actually consume the intermediate
    int_consumed = any(intermediate in set(np.nonzero(rxnMat[r])[0]) 
                    for r in donor_pathway)
    if not int_consumed:
        print(f"ERROR: Donor pathway does not consume the intermediate.", flush=True)
        return None
    
    return donor_pathway

def find_removable_core(net, rxnMat, prodMat, sumRxnVec,
                        nutrientSet, Currency, Core, intermediate):
    '''
    Search for a core metabolite whose production can be removed from the network 
    without interrupting the production of a specified intermediate. 
    '''
    
    removable_core, pruned_net = None, None

    for core in np.random.permutation(Core):
            pruned = remove_core_production(
                net, core, rxnMat, prodMat, sumRxnVec,
                nutrientSet, Currency, Core, [intermediate])
            if pruned is not None:
                removable_core, pruned_net = core, pruned
                break

    if pruned_net is None:
        return None, None
    else:
        return removable_core, pruned_net
    

def build_pathway_pair(net_A, net_B, candidates_A, candidates_B,
                           rxnMat, prodMat, sumRxnVec, nutrientSet, Currency,
                           Core, max_attempts):
    '''
    Searches for a valid cross-feeding pair by repeatedly sampling random
    combinations of candidate intermediates from each network and testing
    feasibility. For each combination, attempts to (1) prune a core metabolite
    from each network without disrupting the chosen intermediate, and (2) build
    a minimal donor pathway from each intermediate to the opposing network's
    missing core. Skips previously attempted intermediate combinations.

    Returns a dict with keys pruned_A, pruned_B, pathway_AB, pathway_BA,
    intermediate_A, intermediate_B, core_A, core_B on success, or None if
    no valid combination is found within max_attempts.
    '''
    attempted = set()
    for attempt in range(1, max_attempts + 1):
        print(f"[Attempt {attempt}/{max_attempts}]", flush=True)
        
        # Try a random combination of intermediates and cores for A and B
        int_A = candidates_A[np.random.randint(len(candidates_A))]
        int_B = candidates_B[np.random.randint(len(candidates_B))]

        # Check if combination has already been attempted
        key = (int(int_A), int(int_B))
        if key in attempted:
            print(f"  Combination already attempted - skipping.", flush=True)
            continue
        attempted.add(key)

        # Attempt to prune core production without interrupting intermediate production
        c_A, p_A = find_removable_core(net_A, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core, int_A)
        if p_A is None:
            print(f"  Cannot prune core {c_A} in net_A without interrupting intermediate {int_A}", flush=True)
            continue

        c_B, p_B = find_removable_core(net_B, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core, int_B)
        if p_B is None:
            print(f"  Cannot prune core {c_B} in net_B without interrupting intermediate {int_B}", flush=True)
            continue

        print(f"  Pairs: A donates {int_A} -> core {c_B} in B | "
              f"B donates {int_B} -> core {c_A} in A", flush=True)

        # Build minimal pathway int_A -> core c_B (using intermediate + currency only)
        path_AB = make_donor_pathway(int_A, c_B, rxnMat, prodMat, sumRxnVec, Currency)
        if path_AB is None or len(path_AB) == 0:
            print(f"  No pathway from {int_A} to core {c_B}.", flush=True)
            continue

        # Build minimal pathway int_B -> core c_A (using intermediate + currency only)
        path_BA = make_donor_pathway(int_B, c_A, rxnMat, prodMat, sumRxnVec, Currency)
        if path_BA is None or len(path_BA) == 0:
            print(f"  No pathway from {int_B} to core {c_A}.", flush=True)
            continue

        if path_AB is not None and path_BA is not None:
            print(f"  Valid pair found! Pathway A --> B has {len(path_AB)} reactions, "
                  f"Pathway B --> A has {len(path_BA)} reactions.", flush=True)

        return {
            'pruned_A': p_A, 
            'pruned_B': p_B,
            'pathway_AB': path_AB, 
            'pathway_BA': path_BA,
            'intermediate_A': int_A, 
            'intermediate_B': int_B,
            'core_A': c_A, 
            'core_B': c_B,}

    print(f"\nERROR: Cross-feeding construction failed after {max_attempts} attempts.",
          flush=True)
    return None


def augment_network(pruned_receiver, donor_pathway, donated_met, rxnMat, prodMat, sumRxnVec, 
                    nutrientSet, Currency, Core, intermediates):
    '''
    Merges a pruned receiver network with a donor pathway, then applies
    alt_randMinNetwork to produce a minimal augmented network that is both
    viable (all core metabolites reachable with donated_met) and obligately
    dependent (not all cores reachable without donated_met). The donated
    metabolite is treated as an additional nutrient only during viability
    checks, not as a base-nutrient supply.

    Returns the pruned augmented reaction array, or None if the batch pruning
    could not reduce the merged network below its initial size (indicating that
    no donor-dependent minimal form exists under the given constraints).
    '''
    combined = np.union1d(pruned_receiver, donor_pathway).astype(int)

    print(f"Augmented receiver size: {len(combined)}", flush=True)

    protected = list(Core)
    if intermediates is not None:
        protected.extend(intermediates)

    rxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    rxnVec[combined] = 1

    augmented_receiver = alt_randMinNetwork(rxnVec, rxnMat, prodMat, sumRxnVec,
                            protected, nutrientSet, Currency, donated_met)
            
    print(f"Pruned modified network: {len(augmented_receiver)}", flush=True)

    if len(augmented_receiver) == len(combined):
        return None
    return augmented_receiver


# ---------------------------------------------------------------------------
# Full obligate cross-feeding pair construction
# ---------------------------------------------------------------------------

def build_crossfeeding_pair(net_A, net_B, rxnMat, prodMat, sumRxnVec,
                             nutrientSet, Currency, Core,
                             use_byproducts, max_attempts, max_runs):
    '''
    Full pipeline for constructing an obligate cross-feeding pair from two
    autonomous networks. Proceeds in four steps: 
    
    (1) Identify candidate secretable intermediates in each network; 
    (2) Search for a valid intermediate/core combination that allows 
    pruning a core from each network and building a donor pathway to 
    it from the partner's intermediate; 
    (3) Augment each pruned receiver with the partner's donor pathway 
    and prune till donor-dependent; 
    (4) Verify that each final network is viable with the donated metabolite 
    but not without it (obligate mutual dependence). If any step fails the procedure 
    retries from step 2 up to max_runs times.

    use_byproducts controls whether candidates are drawn from strict byproducts
    (True) or all non-excluded intermediates (False).

    Returns a dict with keys cross_A, cross_B, A_donated, B_donated suitable
    for passing directly to splitByDemand_crossfeeding, or None on failure.
    '''
    print(f'Network A: {len(net_A)} reactions | Network B: {len(net_B)} reactions', flush=True)

    # Step 1: Identify candidate intermediates from both autonomous networks
    print(f"\n--- Step 1: Identifying candidate intermediates ---", flush=True)
    candidates_A = get_candidates(
        net_A, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        use_byproducts=use_byproducts)

    candidates_B = get_candidates(
        net_B, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        use_byproducts=use_byproducts)

    print(f"Candidates from A: {len(candidates_A)} \nCandidates from B: {len(candidates_B)}", flush=True)

    if len(candidates_A) == 0 or len(candidates_B) == 0:
        print("ERROR: No candidate intermediates found in one or both networks.", flush=True)
        return None
    
    runs = 0

    while True:
        runs += 1
        if runs >= max_runs:
            print(f"\nERROR: Exceeded maximum attempts ({max_runs}) without finding a valid pair.",
                  flush=True)
            return None
        
        # Step 2: Search for valid pathway pairs by pruning on demand
        print(f"\n--- Step 2: Searching for valid pairs (max {max_attempts} attempts) ---",
              flush=True)

        pair = build_pathway_pair(
            net_A, net_B, candidates_A, candidates_B,
            rxnMat, prodMat, sumRxnVec, nutrientSet, Currency,
            Core,
            max_attempts=max_attempts)

        if pair is None:
            return None

        pruned_A      = pair['pruned_A']
        pruned_B      = pair['pruned_B']
        pathway_AB    = pair['pathway_AB']
        pathway_BA    = pair['pathway_BA']
        intermediate_A = pair['intermediate_A']
        intermediate_B = pair['intermediate_B']
        core_A        = pair['core_A']
        core_B        = pair['core_B']

        # Step 3: Augment pruned receivers with donor pathways
        print(f"\n--- Step 3: Augmenting pruned receivers ---", flush=True)

        modified_A = augment_network(
            pruned_A, pathway_BA, intermediate_B,
            rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core,
            [intermediate_A])
        
        modified_B = augment_network(
            pruned_B, pathway_AB, intermediate_A,
            rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core, 
            [intermediate_B])

        if modified_A is None:
            print("ERROR: Could not prune augmented net_A to a donor-dependent state.",
                  flush=True)
            continue
        
        if modified_B is None:
            print("ERROR: Could not prune augmented net_B to a donor-dependent state.",
                  flush=True)
            continue

        # Step 4: Final verification of obligate dependence
        print(f"\n--- Step 4: Final verification of obligate dependence ---", flush=True)

        a_augmented_viable, _ = verify_autonomy(modified_A, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)
        
        b_augmented_viable, _ = verify_autonomy(modified_B, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)
        
        if a_augmented_viable:
            print(f"ERROR: Modified net_A is autonomous without exchange. Retrying from step 2.", flush=True)
            continue
        if b_augmented_viable:
            print(f"ERROR: Modified net_B is autonomous without exchange. Retrying from step 2.", flush=True)
            continue

        a_nutrients = nutrientSet + [intermediate_B]
        a_viable, a_missing = verify_autonomy(
            modified_A, rxnMat, prodMat, sumRxnVec, a_nutrients, Currency, Core)

        b_nutrients = nutrientSet + [intermediate_A]
        b_viable, b_missing = verify_autonomy(
            modified_B, rxnMat, prodMat, sumRxnVec, b_nutrients, Currency, Core)

        if not a_viable:
            print(f"ERROR: Modified net_A not viable with exchange. Missing: {a_missing}",
                  flush=True)
            return None
        if not b_viable:
            print(f"ERROR: Modified net_B not viable with exchange. Missing: {b_missing}",
                  flush=True)
            return None

        print(f"Obligate cross-feeding pair constructed successfully.", flush=True)

        return {
            'auto_A': net_A,
            'auto_B': net_B,
            'cross_A': modified_A,
            'cross_B': modified_B,
            'A_donated': intermediate_A,
            'B_donated': intermediate_B,
            'A_ext_core': core_A,
            'B_ext_core': core_B,
            'pathway_AB': pathway_AB,
            'pathway_BA': pathway_BA,
        }


# ---------------------------------------------------------------------------
# Test Run
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import pickle
    from load_data import *
    import time

    # ── Config ──────────────────────────────────────────────────────────────
    AUTONET_SUBDIR  = "autonets_rs_av1"   # autonets_{source}_av{version}
    AUTONET_FILE    = "P"                  # P (rs is always pruned)
    NET_IDX_A       = 530
    NET_IDX_B       = 825
    # ──────────────────────────────────────────────────────────────────
    from directory_paths import resolve_autonet_path
    with open(resolve_autonet_path(AUTONET_SUBDIR, AUTONET_FILE), "rb") as f:
        all_autonets = pickle.load(f)

    net_A = all_autonets[NET_IDX_A]
    net_B = all_autonets[NET_IDX_B]

    start = time.time()
    crossfeeders = build_crossfeeding_pair(
                             net_A, net_B, rxnMat, prodMat, sumRxnVec,
                             nutrientSet, Currency, Core,
                             use_byproducts=True, max_attempts=10, max_runs=3)
    end = time.time()
    print(f"\nTime taken: {end - start:.2f} seconds")

    if crossfeeders is not None:
        print("\nCross-feeding pair details:")
        print(f"    A needs metabolite {crossfeeders['B_donated']} from B to produce {crossfeeders['A_ext_core']}")
        print(f"    B needs metabolite {crossfeeders['A_donated']} from A to produce {crossfeeders['B_ext_core']}")
        print(f"    A --> B donor pathway has {len(crossfeeders['pathway_AB'])} reactions")
        print(f"    B --> A donor pathway has {len(crossfeeders['pathway_BA'])} reactions")
        print(f"    Net A: {len(crossfeeders['auto_A'])} --> {len(crossfeeders['cross_A'])} reactions")
        print(f"    Net B: {len(crossfeeders['auto_B'])} --> {len(crossfeeders['cross_B'])} reactions")
    else:
        print("\nFailed to construct cross-feeding pair.")