import numpy as np
from satisfiability_check import markSatMetsRxns

def prunedSatsMets(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                   nutrientSet, Currency):
    '''
    Removes reaction(s) remRxns from the current satisfied reaction set satRxns, 
    then re-runs satisfiability propagation from the nutrient+currency
    seed set. Used by isCoreProduced to inspect the removability of reaction(s).

    Accepts a single reaction index or an array of indices as remRxns.

    Returns (satMetVec, satRxnVec) for the pruned network, where
    each vector marks reachable metabolites and reactions respectively.
    '''
    # Creating a temporary set of reactions with some reactions pruned.
    tempSatRxns = np.copy(satRxns)
    tempSatRxns[remRxns] = 0
    
    # Calculating the marked set of reactions from the temporary set.
    return markSatMetsRxns(tempSatRxns, rxnMat, prodMat, sumRxnVec, 
                           nutrientSet, Currency)

#-------------------------------------------------------------------------

def isCoreProduced(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                    nutrientSet, Currency, coreTBP):
    '''
    Tests whether all target core metabolite(s) remain reachable after removing
    reaction(s) remRxns from the network. Used as the acceptance criterion
    throughout every pruning routine: a reaction (or batch) is removable if
    and only if this function returns True.

    coreTBP may be a single metabolite index or a list/array of indices;
    all must remain satisfiable for the function to return True.

    Returns True if all coreTBP metabolites are still produced, False otherwise.
    '''
    tempSatMets, tempSatRxns = prunedSatsMets(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                                              nutrientSet, Currency)

    cores = np.atleast_1d(np.asarray(coreTBP)).ravel()
    return all(tempSatMets[c] for c in cores)