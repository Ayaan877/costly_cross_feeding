import numpy as np
from copy import deepcopy


def isLimiting(tRct, tRxn, m, S, reactants):
    for oRct in reactants:
        if (m[tRct] * S[tRxn, oRct] / S[tRxn, tRct]) > m[oRct]:
            return True


def giveLimitingCurrency(r, tRxn):
    return np.where(r[tRxn] == max(r[tRxn][np.where(r[tRxn] < 0.0)]))[0][0]


def splitByDemand_crossfeeding_alt(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                   rho, pi, nutrientSet, Energy, Currency,
                                   Core, crossPair):
    """
    Coupled two-compartment extension of splitByDemand_alt

    Key behavioural differences vs calculate_autoNet_yield.py:
      - Picks the first NON-limiting reactant (false limiting condition)
      - No nutrient-share persistence across rounds (violates mass conservation)
      - No sumRxnVec == 0 filtering on procRxnVec (fails to filter invalid rxns)
      - Viable if any Core is produced (false viability condition)
    """
    reqKeys = {'cross_A', 'cross_B', 'A_donated', 'B_donated'}
    missing = reqKeys.difference(crossPair.keys())
    if missing:
        raise KeyError(f'crossPair is missing required keys: {", ".join(sorted(missing))}')

    nRxns, nMets = np.shape(stoich_matrix)

    cross_A   = np.array(crossPair['cross_A'], dtype=int)
    cross_B   = np.array(crossPair['cross_B'], dtype=int)
    metA_to_B = int(crossPair['A_donated'])
    metB_to_A = int(crossPair['B_donated'])

    if metA_to_B < 0 or metA_to_B >= nMets:
        raise ValueError('A_donated index is outside stoichiometric bounds.')
    if metB_to_A < 0 or metB_to_A >= nMets:
        raise ValueError('B_donated index is outside stoichiometric bounds.')

    activeA = np.where(cross_A)[0] if cross_A.dtype == bool else cross_A
    activeB = np.where(cross_B)[0] if cross_B.dtype == bool else cross_B
    nA, nB  = len(activeA), len(activeB)

    runningEA, runningBA = 0.0, 0.0
    runningEB, runningBB = 0.0, 0.0
    isCoreProducedA = np.zeros(nMets)
    isCoreProducedB = np.zeros(nMets)
    fluxAtoB, fluxBtoA = 0.0, 0.0

    # Seed both compartments with currency and nutrients.
    metState = np.zeros(2 * nMets)
    seedMets = np.array(list(Currency) + list(nutrientSet), dtype=int)
    if len(seedMets) > 0:
        metState[np.unique(seedMets)]         = 1
        metState[np.unique(seedMets + nMets)] = 1

    # Build coupled matrices: A subnet + B subnet + 2 exchange reactions.
    nPairRxns     = nA + nB + 2
    S             = np.zeros((nPairRxns, 2 * nMets))
    r             = np.zeros((nPairRxns, 2 * nMets))
    rMatPair      = np.zeros((nPairRxns, 2 * nMets))
    pMatPair      = np.zeros((nPairRxns, 2 * nMets))
    sumRxnVecPair = np.zeros(nPairRxns)

    if nA > 0:
        S[:nA, :nMets]        = stoich_matrix[activeA]
        r[:nA, :nMets]        = rho[activeA]
        rMatPair[:nA, :nMets] = rxnMat[activeA]
        pMatPair[:nA, :nMets] = prodMat[activeA]
        sumRxnVecPair[:nA]    = sumRxnVec[activeA]

    if nB > 0:
        S[nA:nA+nB, nMets:]        = stoich_matrix[activeB]
        r[nA:nA+nB, nMets:]        = rho[activeB]
        rMatPair[nA:nA+nB, nMets:] = rxnMat[activeB]
        pMatPair[nA:nA+nB, nMets:] = prodMat[activeB]
        sumRxnVecPair[nA:nA+nB]    = sumRxnVec[activeB]

    rxnAtoB = nA + nB
    rxnBtoA = nA + nB + 1

    S[rxnAtoB, metA_to_B]              = -1
    S[rxnAtoB, nMets + metA_to_B]      =  1
    r[rxnAtoB, metA_to_B]              = -1
    rMatPair[rxnAtoB, metA_to_B]       =  1
    pMatPair[rxnAtoB, nMets + metA_to_B] = 1
    sumRxnVecPair[rxnAtoB]             =  1

    S[rxnBtoA, nMets + metB_to_A]      = -1
    S[rxnBtoA, metB_to_A]              =  1
    r[rxnBtoA, nMets + metB_to_A]      = -1
    rMatPair[rxnBtoA, nMets + metB_to_A] = 1
    pMatPair[rxnBtoA, metB_to_A]       =  1
    sumRxnVecPair[rxnBtoA]             =  1

    procRxnVec = ((np.dot(rMatPair, metState != 0) - sumRxnVecPair) == 0) * 1
    # Note: no sumRxnVecPair == 0 filtering

    isChecked = np.zeros(nPairRxns)

    mask = np.abs(np.sum(r, axis=0)) != 0
    shareMatrix = np.zeros(S.shape, dtype=float)
    shareMatrix[:, np.where(mask)[0]] = ((r * metState)[:, mask] /
                                         np.abs(np.sum(r, axis=0))[mask])

    currencyAB = list(Currency) + [c + nMets for c in Currency]
    if currencyAB:
        shareMatrix[:, currencyAB] = -1

    trackedEnergy = set(Energy)
    trackedCore   = set(Core)
    trackedMetsA  = np.array(list(trackedCore | trackedEnergy), dtype=int)
    trackedMetsB  = trackedMetsA + nMets
    trackedMets   = np.append(trackedMetsA, trackedMetsB)

    while procRxnVec.any():
        prodState = np.zeros(2 * nMets)

        for thisRxn in np.where(procRxnVec)[0]:
            allowedRct = []
            isChecked[thisRxn] = 1

            rs        = np.where(rMatPair[thisRxn])[0]
            ps        = np.where(pMatPair[thisRxn])[0]
            reactants = [tR for tR in rs if tR not in currencyAB]
            products  = [tP for tP in ps if tP not in currencyAB]

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

            if thisRxn == rxnAtoB:
                fluxAtoB += -shareMatrix[thisRxn, limRct]
            elif thisRxn == rxnBtoA:
                fluxBtoA += -shareMatrix[thisRxn, limRct]

            mets = np.append(rs, ps)
            for thisMet in mets[np.where(np.isin(mets, trackedMets))]:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]

                if thisMet in trackedEnergy:
                    runningEA += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in [e + nMets for e in trackedEnergy]:
                    runningEB += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in trackedCore:
                    if thisMet in ps:
                        isCoreProducedA[Core] = 1
                    runningBA += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in [c + nMets for c in trackedCore]:
                    if thisMet in ps:
                        isCoreProducedB[Core] = 1
                    runningBB += shareMatrix[thisRxn, limRct] * ratio

            for thisMet in reactants:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                shareMatrix[thisRxn, thisMet] -= shareMatrix[thisRxn, limRct] * ratio

        # Redistribution: update ALL columns with remaining demand
        r[np.where(isChecked)] = 0
        mask = np.abs(np.sum(r, axis=0)) != 0
        shareMatrix[:, np.where(mask)[0]] = ((r * prodState)[:, mask] /
                                             np.abs(np.sum(r, axis=0))[mask])
        if currencyAB:
            shareMatrix[:, currencyAB] = -1
        # Note: no nutrient-share persistence across rounds

        procRxnVec = ((np.dot(rMatPair, np.sum(shareMatrix, axis=0) != 0) -
                       sumRxnVecPair) == 0) * 1
        procRxnVec[np.where(isChecked)] = 0
        # Note: no sumRxnVecPair == 0 filtering

    statusA = bool(isCoreProducedA[Core].all())
    statusB = bool(isCoreProducedB[Core].all())

    return {
        'E_A':         float(runningEA),
        'B_A':         float(runningBA),
        'viable_A':    statusA,
        'E_B':         float(runningEB),
        'B_B':         float(runningBB),
        'viable_B':    statusB,
        'pair_viable': bool(statusA and statusB),
        'flux_A_to_B': float(fluxAtoB),
        'flux_B_to_A': float(fluxBtoA),
    }
