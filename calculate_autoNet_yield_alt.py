import numpy as np
from copy import deepcopy


def isLimiting(tRct, tRxn, m, S, reactants):
    for oRct in reactants:
        if (m[tRct] * S[tRxn, oRct] / S[tRxn, tRct]) > m[oRct]:
            return True


def giveLimitingCurrency(r, tRxn):
    return np.where(r[tRxn] == max(r[tRxn][np.where(r[tRxn] < 0.0)]))[0][0]


def splitByDemand_alt(stoich_matrix, rxnMat, prodMat, sumRxnVec, rho, pi,
                      nutrientSet, Energy, Currency, Core, orgRxns):
    """
    Key behavioural differences vs calculate_autoNet_yield.py:
      - Picks the first NON-limiting reactant (false limiting condition)
      - No nutrient-share persistence across rounds (violates mass conservation)
      - No sumRxnVec == 0 filtering on procRxnVec (fails to filter invalid rxns)
      - Viable if any Core is produced (false viability condition)
      - Returns (-1.0, -1.0, False) for non-viable networks (returns float for non-viable nets)
    """
    runningE, runningB = 0.0, 0.0

    nMets      = len(np.transpose(stoich_matrix))
    activeRxns = np.array(orgRxns, dtype=int)
    nActive    = len(activeRxns)
    isCoreProduced = np.zeros(nMets)

    metState = np.zeros(nMets)
    metState[Currency + nutrientSet] = 1

    r    = np.copy(rho[activeRxns])
    S    = np.copy(stoich_matrix[activeRxns])
    rMat = np.copy(rxnMat[activeRxns])
    pMat = np.copy(prodMat[activeRxns])
    sumRxnVecActive = sumRxnVec[activeRxns]

    procRxnVec = ((np.dot(rMat, metState != 0) - sumRxnVecActive) == 0) * 1

    isChecked = np.zeros(nActive)

    mask = np.abs(np.sum(r, axis=0)) != 0
    shareMatrix = np.zeros(S.shape, dtype=float)
    shareMatrix[:, np.where(mask)[0]] = ((r * metState)[:, mask] /
                                         np.abs(np.sum(r, axis=0))[mask])
    shareMatrix[:, Currency] = -1

    while procRxnVec.any():
        prodState = np.zeros(nMets)

        for thisRxn in np.where(procRxnVec)[0]:
            allowedRct = []
            isChecked[thisRxn] = 1

            rs        = np.where(rMat[thisRxn])[0]
            ps        = np.where(pMat[thisRxn])[0]
            reactants = [tR for tR in rs if tR not in Currency]
            products  = [tP for tP in ps if tP not in Currency]

            # Pick the first non-limiting reactant
            for thisReactant in reactants:
                if not isLimiting(thisReactant, thisRxn,
                                   shareMatrix[thisRxn], S, reactants):
                    allowedRct.append(thisReactant)
                    limRct = deepcopy(thisReactant)
                    break

            if not allowedRct:
                limRct = giveLimitingCurrency(r, thisRxn)

            for thisMet in products:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                prodState[thisMet] += shareMatrix[thisRxn, limRct] * ratio

            mets = np.append(rs, ps)
            for thisMet in mets[np.where(np.isin(mets, np.array(Core + Energy)))]:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]

                if thisMet in Energy:
                    runningE += shareMatrix[thisRxn, limRct] * ratio

                # Viable if any Core is produced, even if not all are
                elif thisMet in Core:
                    if thisMet in ps:
                        isCoreProduced[Core] = 1
                    runningB += shareMatrix[thisRxn, limRct] * ratio

            for thisMet in reactants:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                shareMatrix[thisRxn, thisMet] -= shareMatrix[thisRxn, limRct] * ratio
        
        # Redistribution: update ALL columns with remaining demand
        r[np.where(isChecked)] = 0
        mask = np.abs(np.sum(r, axis=0)) != 0
        shareMatrix[:, np.where(mask)[0]] = ((r * prodState)[:, mask] /
                                             np.abs(np.sum(r, axis=0))[mask])
        shareMatrix[:, Currency] = -1
        # Note: no nutrient-share persistence across rounds

        procRxnVec = ((np.dot(rMat, np.sum(shareMatrix, axis=0) != 0) -
                       sumRxnVecActive) == 0) * 1
        procRxnVec[np.where(isChecked)] = 0
        # Note: no sumRxnVec == 0 filtering

    if isCoreProduced[Core].all():
        return runningE, runningB, True
    else:
        return -1.0, -1.0, False
