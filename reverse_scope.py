import numpy as np
from satisfiability_check import markSatMetsRxns, make_sparse

def giveRevScope(rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, coreTBP):
    '''
    Discovers the complete non-minimal set of reactions that could contribute 
    to producing a target core metabolite from the given nutrient and currency 
    seed set via reverse scope expansion. Starting from coreTBP, the algorithm
    iteratively marks reactions that produce "frontier" metabolites and their
    required reactants, expanding backward until the nutrient set fully
    satisfies all marked reactions. Satisfiability is verified at each step
    with markSatMetsRxns (forward scope) to ensure that only reactions genuinely 
    reachable from the seeds are included.

    The returned subgraph is non-minimal; it includes all reactions in the
    reverse scope, many of which may be redundant. It should be passed to
    a pruning routine (randMinNetwork or similar) to obtain a minimal pathway.

    coreTBP may be a single metabolite index or an array of indices.

    Returns (satMets, satRxns), binary vectors of length nMets and nRxns
    marking the reachable metabolites and reactions in the reverse scope.
    Raises ValueError if the target metabolite cannot be reached from the
    nutrient set at all.
    '''

    # Convert to sparse once if needed (cached for repeated calls).
    sp_rxnMat = make_sparse(rxnMat)
    sp_prodMat = make_sparse(prodMat)
    sp_rxnMat_T = sp_rxnMat.T.tocsr()

    n_rxns, n_mets = sp_rxnMat.get_shape()

    # Initializing all the vectors to propagate the satisfied subgraph search.
    seedVec = np.zeros(n_mets)
    rxnProc = np.zeros(n_rxns)
    seedVec[coreTBP] = 1
    currScopeMets = np.copy(seedVec)
    prevScopeMets = np.copy(seedVec)
    deltaMetVec = np.copy(seedVec)

    while True:
        prevScopeMets = np.logical_or(currScopeMets, prevScopeMets)

        # Propagating reverse scope.
        prod_delta = np.asarray(sp_prodMat.dot(deltaMetVec)).ravel()
        rxnProc = (prod_delta + rxnProc > 0) * 1
        currScopeMets = (np.asarray(sp_rxnMat_T.dot(prod_delta)).ravel() > 0) * 1
        currScopeMets = np.logical_xor(currScopeMets, np.logical_and(currScopeMets, prevScopeMets))

        # Marking the satisfied metabolites and reactions.
        satMets, satRxns = markSatMetsRxns(rxnProc, sp_rxnMat, sp_prodMat, sumRxnVec, nutrientSet, Currency)

        # If core has been reached, checking if everything is marked, then returning.
        if np.array_equal(satRxns, rxnProc):
            return satMets, satRxns

        # Calculating the new metabolites that need to be produced
        deltaMetVec = np.logical_xor(currScopeMets, 
                                     np.logical_and(currScopeMets, satMets)) * 1

        # If the full reverse scope has been reached, stopping and returning the marked set.
        if set(np.nonzero(currScopeMets)[0]).issubset(set(np.nonzero(prevScopeMets)[0])):

            # Verify all requested cores are satisfied.
            cores = np.atleast_1d(np.asarray(coreTBP)).ravel()
            unsatisfied = [int(c) for c in cores if not satMets[c]]
            if unsatisfied:
                raise ValueError(f"Reverse scope failed to satisfy target(s): {unsatisfied}")

            return satMets, satRxns      