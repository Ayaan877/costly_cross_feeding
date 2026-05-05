import numpy as np
from satisfiability_check import markSatMetsRxns

def verify_autonomy(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core):
    '''
    Checks whether a network can produce all core biomass precursors starting
    from the given nutrient set and currency metabolites, using reachability
    (satisfiability) alone. 
    
    Note: This is a necessary but not sufficient condition for viability. 
    This only checks if all core metabolites are reachable in principle, 
    but their flux may still be limited by stoichiometric constraints.

    Returns (is_autonomous, missing_cores) where is_autonomous is True if all
    Core metabolites are reachable and missing_cores lists any that are not.
    '''
    rxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    rxnVec[net] = 1
    satMets, _ = markSatMetsRxns(rxnVec, rxnMat, prodMat,
                                 sumRxnVec, nutrientSet, Currency)

    missing_cores = [c for c in Core if not satMets[c]]

    return len(missing_cores) == 0, missing_cores