import numpy as np
from copy import deepcopy

def isLimiting(tRct, tRxn, m, S, reactants):
    for oRct in reactants:
        if (m[tRct] * S[tRxn, oRct] / S[tRxn,tRct]) > m[oRct]:
            return True

def giveLimitingCurrency(r, tRxn):
    return np.where(r[tRxn] == max(r[tRxn][np.where(r[tRxn] < 0.0)]))[0][0]

def splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec, rho, pi,
                 nutrientSet, Energy, Currency, Core, orgRxns, nutrientSupply=None):
    '''
    Takes a set of reactions that form a complete network (orgRxns) and calculates the energy 
    and biomass yields based on stoichiometric demand. Flux is allocated according to the limiting 
    reactant for each reaction, and the total demand across all reactions, including downstream
    demand. Flux accumulates across rounds to ensure mass conservation.
    '''
    # Initializing yield counters for E and B.
    runningE, runningB = 0.0, 0.0

    # Getting the compact set of active reactions in this network.
    nMets = len(np.transpose(stoich_matrix))
    activeRxns = np.array(orgRxns, dtype = int)
    nActive = len(activeRxns)
    isCoreProduced = np.zeros(nMets)

    # Constructing a new metabolite state vector.
    metState = np.zeros(nMets)
    metState[ Currency + nutrientSet ] = 1
    if nutrientSupply:
        for met, supply in nutrientSupply.items():
            metState[met] = supply

    # Getting compact copies of the relevant matrices (active reactions only).
    r = np.copy(rho[activeRxns])
    S = np.copy(stoich_matrix[activeRxns])
    rMat = np.copy(rxnMat[activeRxns])
    pMat = np.copy(prodMat[activeRxns])
    sumRxnVecActive = sumRxnVec[activeRxns]

    # Figuring out which reactions can be performed at this step.
    procRxnVec = ((np.dot(rMat, metState != 0) - sumRxnVecActive) == 0) * 1
    procRxnVec[sumRxnVecActive == 0] = 0

    # Continuing calculation till no more reactions can be performed
    isChecked = np.zeros(nActive)

    # Computing which reaction gets what share of which reactants.
    mask = np.abs(np.sum(r, axis = 0)) != 0
    shareMatrix = np.zeros(S.shape, dtype=float)
    shareMatrix[:, np.where(mask)[0]] = ((r * metState)[:, mask] / 
                                         np.abs(np.sum(r, axis = 0))[mask])
    shareMatrix[:, Currency] = -1

    # Saving initial total demand for nutrients to maintain their share across rounds.
    totalInitialDemand = np.abs(np.sum(r, axis = 0))

    while procRxnVec.any():

        # Initializing the product metabolite state vector.
        prodState = np.zeros(len(np.transpose(stoich_matrix)))

        # Updating states after all accomplishable reactions.
        for thisRxn in np.where(procRxnVec)[0]:
            # Checking if found a usable reactant.
            allowedRct = []
            isChecked[thisRxn] = 1

            # Getting the reactants and products of this reaction, except currency.
            rs, ps = np.where(rMat[thisRxn])[0], np.where(pMat[thisRxn])[0]
            reactants = [tR for tR in rs if tR not in Currency]
            products = [tP for tP in ps if tP not in Currency]

            # Checking for limiting reactants.
            for thisReactant in reactants:
                if isLimiting(thisReactant, thisRxn, shareMatrix[thisRxn], S, reactants):
                    allowedRct.append(thisReactant)
                    limRct = deepcopy(thisReactant)
                    break

            # If nothing is limiting, everything gets used.
            if not allowedRct:
                if reactants:
                    limRct = reactants[0]
                else:
                    limRct = giveLimitingCurrency(r, thisRxn)

            # Updating metabolite amounts post reaction.
            for thisMet in products:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                prodState[thisMet] += shareMatrix[thisRxn, limRct] * ratio

            mets = np.append(rs, ps)
            for thisMet in mets[np.where(np.isin(mets, np.array(Core + Energy)))]:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]

                # Updating E and B if these metabolites are produced or consumed.
                if thisMet in Energy:
                    runningE += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in Core:
                    if thisMet in ps:
                        isCoreProduced[thisMet] = 1
                    runningB += shareMatrix[thisRxn, limRct] * ratio

            # Updating metabolite amounts post reaction.
            for thisMet in reactants:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                shareMatrix[thisRxn, thisMet] -= shareMatrix[thisRxn, limRct] * ratio

        # Redistributing the produced metabolites among reactions by demand.
        r[np.where(isChecked)] = 0
        totalDemand = np.abs(np.sum(r, axis = 0))
        newProdCols = np.where((prodState != 0) & (totalDemand != 0))[0]
        if len(newProdCols) > 0:
            shareMatrix[:, newProdCols] += ((r * prodState)[:, newProdCols] /
                                            totalDemand[newProdCols])
        shareMatrix[:, Currency] = -1

        # Maintaining nutrient shares across rounds using initial total demand.
        for met in nutrientSet:
            if totalInitialDemand[met] > 0:
                supply = nutrientSupply.get(met, 1.0) if nutrientSupply else 1.0
                shareMatrix[:, met] = supply * r[:, met] / totalInitialDemand[met]

        # Recalculating performable reactions.
        procRxnVec = ((np.dot(rMat, np.sum(shareMatrix, axis = 0) != 0) - sumRxnVecActive) == 0) * 1
        procRxnVec[np.where(isChecked)] = 0
        procRxnVec[sumRxnVecActive == 0] = 0

    if isCoreProduced[Core].all():
        return runningE, runningB, True
    else:
        return None, None, False
