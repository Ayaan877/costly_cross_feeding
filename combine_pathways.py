import numpy as np
from prune_check import isCoreProduced

def buildAutonomousNetwork(pathway_list, rxnMat, prodMat, sumRxnVec, 
                           nutrientSet, Currency, coreTBPs, prune, rng=None):
    '''
    Constructs an autonomous network capable of producing all 8 core biomass
    precursors by taking the union of 8 individual minimal pathways (one per
    core target). If prune=True, the merged union is further reduced by a
    single-reaction sweep that removes any reaction whose deletion does not
    prevent production of any core metabolite, producing a minimal autonomous
    network. If prune=False, the raw union is returned.

    pathway_list is a list of 8 arrays of reaction indices, one pathway per
    core target in coreTBPs. coreTBPs is an array of 8 metabolite indices.

    Returns a numpy array of reaction indices for the (optionally minimal)
    autonomous network.
    '''

    if rng is None:
        rng = np.random.default_rng()

    # Union of reactions from 8 pathways
    combined_rxns = set()
    for pathway in pathway_list:
        combined_rxns.update(pathway)

    combined_rxns = np.array(list(combined_rxns), dtype=int)

    satRxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    satRxnVec[combined_rxns] = 1
    
    if prune:
        currSatRxnVec = np.copy(satRxnVec)

        while True:
            currSatRxns = np.nonzero(currSatRxnVec)[0]
            removed_any = False

            for rxn in rng.permutation(currSatRxns):
                if isCoreProduced(rxn, currSatRxnVec, rxnMat, prodMat,
                                    sumRxnVec, nutrientSet, Currency, coreTBPs):
                    currSatRxnVec[rxn] = 0
                    removed_any = True

            if not removed_any:
                break

        return np.nonzero(currSatRxnVec)[0]
    else:
        return combined_rxns