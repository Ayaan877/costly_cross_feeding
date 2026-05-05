import numpy as np
from prune_check import isCoreProduced

def randMinNetwork(satRxnVec, rxnMat, prodMat, sumRxnVec,
                   Core, nutrientSet, Currency, rng=None, init_frac=0.5):
    '''
    Prunes a set of reactions to a minimal network using an adaptive
    two-phase strategy. This hybrid approach is substantially faster than 
    pure single-reaction sweeping on large initial networks.
    
    Batch phase: Iteratively attempts to remove random batches of reactions, 
    accepting any removal that preserves production of all core metabolites.
    Failed attempts adaptively shrink the batch fraction until batch <= 1. 
    
    Single phase: Once batch size reaches 1, reactions are randomly permuted 
    and tested for individual removability until no more can be removed.
    
    Returns a numpy array of reaction indices for the minimal network.
    '''
    # print(f"Batch pruning network...", flush=True)
    if rng is None:
        rng = np.random.default_rng()

    currSatRxnVec = np.copy(satRxnVec)
    fail_count = 0
    max_fails = 3

    batch_frac = init_frac

    while True:
        currSatRxns = np.nonzero(currSatRxnVec)[0]
        n_curr = len(currSatRxns)

        if n_curr == 0:
            break

        batch_size = max(1, int(n_curr * batch_frac))

        # ------------------------------------------------
        # Batchwise Reaction Removal Phase
        # ------------------------------------------------
        if batch_size > 1:
            batch = rng.choice(currSatRxns, size=batch_size, replace=False)

            if isCoreProduced(batch, currSatRxnVec, rxnMat, prodMat, 
                              sumRxnVec, nutrientSet, Currency, Core):

                currSatRxnVec[batch] = 0
                fail_count = 0

            else:
                fail_count += 1
                if fail_count >= max_fails:
                    batch_frac /= 1.1
                    fail_count = 0

        # ------------------------------------------------
        # Single Reaction Systematic Removal Phase
        # ------------------------------------------------
        else:
            removed_any = False

            for rxn in rng.permutation(currSatRxns):
                if isCoreProduced([rxn], currSatRxnVec, rxnMat, prodMat,
                                  sumRxnVec, nutrientSet, Currency, Core):

                    currSatRxnVec[rxn] = 0
                    removed_any = True

            if not removed_any:
                # print(f'Minimal network size = {len(currSatRxns)}', flush=True)
                break
    
    return np.nonzero(currSatRxnVec)[0]



def alt_randMinNetwork(satRxnVec, rxnMat, prodMat, sumRxnVec,
                   protected_mets, nutrientSet, Currency, donated_met, 
                   rng=None, init_frac=0.5):
    '''
    Variant of randMinNetwork that enforces obligate dependence on a donated
    metabolite throughout pruning. Uses the same adaptive batch-then-single
    strategy, but a reaction removal is accepted only when two conditions hold
    simultaneously: 
    (1) All protected metabolites are producible with the donated metabolite 
    present in the nutrient set (viability),
    (2) Not all protected metabolites are producible in the absence of the 
    donated metabolite (dependency).

    protected_mets lists the metabolite indices that must remain reachable
    (typically all Core metabolites plus the organism's own secreted intermediate).
    donated_met is the externally supplied metabolite that must be required 
    for viability (typically the partner's secreted intermediate).

    Returns a numpy array of reaction indices for the minimal, donor-dependent network.
    '''
    print(f"Batch pruning network while preserving dependence on intermediate {donated_met}...", flush=True)
    if rng is None:
        rng = np.random.default_rng()

    augmented_nutrients = list(nutrientSet) + [donated_met]
    base_nutrients = nutrientSet

    currSatRxnVec = np.copy(satRxnVec)
    fail_count = 0
    max_fails = 3

    batch_frac = init_frac

    while True:
        currSatRxns = np.nonzero(currSatRxnVec)[0]
        n_curr = len(currSatRxns)

        if n_curr == 0:
            break

        batch_size = max(1, int(n_curr * batch_frac))

        # ------------------------------------------------
        # Batchwise Reaction Removal Phase
        # ------------------------------------------------
        if batch_size > 1:
            batch = rng.choice(currSatRxns, size=batch_size, replace=False)

            viable   = isCoreProduced(batch, currSatRxnVec, rxnMat, prodMat,
                                      sumRxnVec, augmented_nutrients, Currency, protected_mets)
            dependent = not isCoreProduced(batch, currSatRxnVec, rxnMat, prodMat,
                                           sumRxnVec, base_nutrients, Currency, protected_mets)

            if viable and dependent:
                currSatRxnVec[batch] = 0
                fail_count = 0
            else:
                fail_count += 1
                if fail_count >= max_fails:
                    batch_frac /= 1.05
                    fail_count = 0

        # ------------------------------------------------
        # Single Reaction Systematic Removal Phase
        # ------------------------------------------------
        else:
            removed_any = False

            for rxn in rng.permutation(currSatRxns):
                viable   = isCoreProduced([rxn], currSatRxnVec, rxnMat, prodMat,
                                          sumRxnVec, augmented_nutrients, Currency, protected_mets)
                dependent = not isCoreProduced([rxn], currSatRxnVec, rxnMat, prodMat,
                                               sumRxnVec, base_nutrients, Currency, protected_mets)

                if viable and dependent:
                    currSatRxnVec[rxn] = 0
                    removed_any = True

            if not removed_any:
                # print(f'Minimal network size = {len(currSatRxns)}', flush=True)
                break
    
    return np.nonzero(currSatRxnVec)[0]