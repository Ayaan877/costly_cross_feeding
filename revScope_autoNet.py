from reverse_scope import giveRevScope
from batch_pruning import randMinNetwork

def revScopeAutoNet(rxnMat, prodMat, sumRxnVec,
                            nutrientSet, Currency, coreTBPs):
    '''
    Generates a minimal autonomous network that produces all target core
    metabolites from the nutrient and currency seed set. Runs reverse scope
    expansion simultaneously over all coreTBPs to obtain the initial non-minimal
    subgraph, then applies batch pruning to reduce it to a minimal network.
    
    coreTBPs is a list of metabolite indices for all core biomass precursors.

    Returns a numpy array of reaction indices for the minimal autonomous network.
    '''

    satMets, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec,
                                    nutrientSet, Currency, coreTBPs)
    return randMinNetwork(satRxns, rxnMat, prodMat, sumRxnVec,
                          coreTBPs, nutrientSet, Currency)


# ---------------------------------------------------------------
# Test Run
# ---------------------------------------------------------------

if __name__ == "__main__":
    from load_data import *
    import numpy as np
    import time
    from satisfiability_check import markSatMetsRxns

    start = time.time()
    net = revScopeAutoNet(rxnMat, prodMat, sumRxnVec,
                                  nutrientSet, Currency, Core)
    satRxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    satRxnVec[net] = 1
    satMets, satRxns = markSatMetsRxns(satRxnVec, rxnMat, prodMat, sumRxnVec,
                             nutrientSet, Currency)
    
    viable = all(satMets[c] for c in Core)
    
    print(f"Generated autonet with {len(net)} reactions: {net}")
    print(f"Network viability: {viable}")
    print(f"Time taken: {(time.time() - start)} seconds")
