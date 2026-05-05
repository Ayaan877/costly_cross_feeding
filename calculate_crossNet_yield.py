import numpy as np
from copy import deepcopy

def isLimiting(tRct, tRxn, m, S, reactants):
    for oRct in reactants:
        if (m[tRct] * S[tRxn, oRct] / S[tRxn, tRct]) > m[oRct]:
            return True

def giveLimitingCurrency(r, tRxn):
    return np.where(r[tRxn] == max(r[tRxn][np.where(r[tRxn] < 0.0)]))[0][0]

def splitByDemand_crossfeeding(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                  rho, pi, nutrientSet, Energy, Currency,
                                  Core, crossPair):
    '''
    Calculates individual energy and biomass yields for a cross-feeding pair
    by solving both metabolic networks simultaneously in two coupled compartments.
    Builds a doubled metabolite space (2 x nMets indices: A in [0, nMets), B in
    [nMets, 2 x nMets)) and adds two explicit exchange reactions that transfer
    A_donated from compartment A to B and B_donated from B to A. Flux is then
    allocated round-by-round using the same stoichiometric demand-splitting
    approach as splitByDemand: metabolite supply is split proportional to
    stoichiometric demand across all consuming reactions, regardless of whether
    they are currently available to fire. Each round, all reactions whose
    allocated reactant shares are sufficient execute; currency metabolites are
    treated as unlimited.

    crossPair is a dict with keys cross_A (reaction indices for organism A),
    cross_B (reaction indices for organism B), A_donated (metabolite index
    secreted by A), and B_donated (metabolite index secreted by B). Exchange
    metabolites start at zero supply, so mutually deadlocked pairs (neither
    organism can fire before receiving the partner's metabolite) will produce
    no exchange flux and be returned as non-viable.

    Returns a dict with keys E_A, B_A, viable_A, E_B, B_B, viable_B,
    pair_viable, flux_A_to_B, flux_B_to_A. Yields are per unit primary nutrient.
    '''
    reqKeys = {'cross_A', 'cross_B', 'A_donated', 'B_donated'}
    missing = reqKeys.difference(crossPair.keys())
    if missing:
        missingStr = ', '.join(sorted(missing))
        raise KeyError(f'crossPair is missing required keys: {missingStr}')

    nRxns, nMets = np.shape(stoich_matrix)

    cross_A = np.array(crossPair['cross_A'], dtype = int)
    cross_B = np.array(crossPair['cross_B'], dtype = int)
    metA_to_B = int(crossPair['A_donated'])
    metB_to_A = int(crossPair['B_donated'])

    if metA_to_B < 0 or metA_to_B >= nMets:
        raise ValueError('A_donated index is outside stoichiometric bounds.')
    if metB_to_A < 0 or metB_to_A >= nMets:
        raise ValueError('B_donated index is outside stoichiometric bounds.')

    activeA = np.where(cross_A)[0] if cross_A.dtype == bool else cross_A
    activeB = np.where(cross_B)[0] if cross_B.dtype == bool else cross_B
    nA = len(activeA)
    nB = len(activeB)

    # Initializing yield counters and production flags for each network.
    runningEA, runningBA = 0.0, 0.0
    runningEB, runningBB = 0.0, 0.0
    isCoreProducedA = np.zeros(nMets)
    isCoreProducedB = np.zeros(nMets)

    fluxAtoB, fluxBtoA = 0.0, 0.0

    # Constructing metabolite state for coupled A/B compartments.
    metState = np.zeros(2 * nMets)
    seedMets = np.array(list(Currency) + list(nutrientSet), dtype = int)
    if len(seedMets) > 0:
        metState[np.unique(seedMets)] = 1
        metState[np.unique(seedMets + nMets)] = 1

    # Building coupled matrices: A subnet + B subnet + 2 exchange links.
    nPairRxns = nA + nB + 2
    S = np.zeros((nPairRxns, 2 * nMets))
    r = np.zeros((nPairRxns, 2 * nMets))
    rMatPair = np.zeros((nPairRxns, 2 * nMets))
    pMatPair = np.zeros((nPairRxns, 2 * nMets))
    sumRxnVecPair = np.zeros(nPairRxns)

    if nA > 0:
        S[:nA, :nMets] = stoich_matrix[activeA]
        r[:nA, :nMets] = rho[activeA]
        rMatPair[:nA, :nMets] = rxnMat[activeA]
        pMatPair[:nA, :nMets] = prodMat[activeA]
        sumRxnVecPair[:nA] = sumRxnVec[activeA]

    if nB > 0:
        S[nA:nA + nB, nMets:] = stoich_matrix[activeB]
        r[nA:nA + nB, nMets:] = rho[activeB]
        rMatPair[nA:nA + nB, nMets:] = rxnMat[activeB]
        pMatPair[nA:nA + nB, nMets:] = prodMat[activeB]
        sumRxnVecPair[nA:nA + nB] = sumRxnVec[activeB]

    # Adding exchange reaction rows.
    rxnAtoB = nA + nB
    rxnBtoA = nA + nB + 1

    S[rxnAtoB, metA_to_B] = -1
    S[rxnAtoB, nMets + metA_to_B] = 1
    r[rxnAtoB, metA_to_B] = -1
    rMatPair[rxnAtoB, metA_to_B] = 1
    pMatPair[rxnAtoB, nMets + metA_to_B] = 1
    sumRxnVecPair[rxnAtoB] = 1

    S[rxnBtoA, nMets + metB_to_A] = -1
    S[rxnBtoA, metB_to_A] = 1
    r[rxnBtoA, nMets + metB_to_A] = -1
    rMatPair[rxnBtoA, nMets + metB_to_A] = 1
    pMatPair[rxnBtoA, metB_to_A] = 1
    sumRxnVecPair[rxnBtoA] = 1

    # Figuring out which reactions can be performed at this step.
    procRxnVec = ((np.dot(rMatPair, metState != 0) - sumRxnVecPair) == 0) * 1
    procRxnVec[sumRxnVecPair == 0] = 0

    # Continuing calculation till no more reactions can be performed.
    isChecked = np.zeros(nPairRxns)

    # Computing which reaction gets what share of which reactants.
    mask = np.abs(np.sum(r, axis = 0)) != 0
    shareMatrix = np.zeros(S.shape, dtype=float)
    shareMatrix[:, np.where(mask)[0]] = ((r * metState)[:, mask] /
                                         np.abs(np.sum(r, axis = 0))[mask])

    currencyAB = list(Currency) + [curr + nMets for curr in Currency]
    if len(currencyAB) > 0:
        shareMatrix[:, currencyAB] = -1

    # Saving initial total demand for nutrients to maintain their share across rounds.
    totalInitialDemandPair = np.abs(np.sum(r, axis = 0))
    nutrientColsAB = list(nutrientSet) + [m + nMets for m in nutrientSet]
    trackedEnergy = set(Energy)
    trackedCore = set(Core)
    trackedMetsA = np.array(list(trackedCore.union(trackedEnergy)), dtype = int)
    trackedMetsB = trackedMetsA + nMets
    trackedMets = np.append(trackedMetsA, trackedMetsB)

    while procRxnVec.any():

        # Initializing the product metabolite state vector.
        prodState = np.zeros(2 * nMets)

        # Updating states after all accomplishable reactions.
        for thisRxn in np.where(procRxnVec)[0]:
            # Checking if found a usable reactant.
            allowedRct = []
            isChecked[thisRxn] = 1

            # Getting reactants and products of this reaction, except currency.
            rs, ps = np.where(rMatPair[thisRxn])[0], np.where(pMatPair[thisRxn])[0]
            reactants = [tR for tR in rs if tR not in currencyAB]
            products = [tP for tP in ps if tP not in currencyAB]

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

            if thisRxn == rxnAtoB:
                fluxAtoB += -shareMatrix[thisRxn, limRct]
            elif thisRxn == rxnBtoA:
                fluxBtoA += -shareMatrix[thisRxn, limRct]

            mets = np.append(rs, ps)
            for thisMet in mets[np.where(np.isin(mets, trackedMets))]:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]

                # Updating E and B for network A or B.
                if thisMet in trackedEnergy:
                    runningEA += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in [e + nMets for e in trackedEnergy]:
                    runningEB += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in trackedCore:
                    if thisMet in ps:
                        isCoreProducedA[thisMet] = 1
                    runningBA += shareMatrix[thisRxn, limRct] * ratio
                elif thisMet in [c + nMets for c in trackedCore]:
                    if thisMet in ps:
                        isCoreProducedB[thisMet - nMets] = 1
                    runningBB += shareMatrix[thisRxn, limRct] * ratio

            # Updating metabolite amounts post reaction.
            for thisMet in reactants:
                ratio = S[thisRxn, thisMet] / S[thisRxn, limRct]
                shareMatrix[thisRxn, thisMet] -= shareMatrix[thisRxn, limRct] * ratio

        # Redistributing produced metabolites among reactions by demand.
        r[np.where(isChecked)] = 0
        totalDemandPair = np.abs(np.sum(r, axis = 0))
        newProdCols = np.where((prodState != 0) & (totalDemandPair != 0))[0]
        if len(newProdCols) > 0:
            shareMatrix[:, newProdCols] += ((r * prodState)[:, newProdCols] /
                                            totalDemandPair[newProdCols])
        if len(currencyAB) > 0:
            shareMatrix[:, currencyAB] = -1

        # Maintaining nutrient shares across rounds using initial total demand.
        for col in nutrientColsAB:
            if totalInitialDemandPair[col] > 0:
                shareMatrix[:, col] = r[:, col] / totalInitialDemandPair[col]

        # Recalculating performable reactions.
        procRxnVec = ((np.dot(rMatPair, np.sum(shareMatrix, axis = 0) != 0) -
                       sumRxnVecPair) == 0) * 1
        procRxnVec[np.where(isChecked)] = 0
        procRxnVec[sumRxnVecPair == 0] = 0

    statusA = bool(isCoreProducedA[Core].all())
    statusB = bool(isCoreProducedB[Core].all())

    return {
        'E_A': float(runningEA),
        'B_A': float(runningBA),
        'viable_A': statusA,
        'E_B': float(runningEB),
        'B_B': float(runningBB),
        'viable_B': statusB,
        'pair_viable': bool(statusA and statusB),
        'flux_A_to_B': float(fluxAtoB),
        'flux_B_to_A': float(fluxBtoA),
    }


if __name__ == "__main__":
    import time
    from load_data import *
    from load_networks import load_minpaths
    from crossfeeding_minPaths import build_crossfeeding_pair_from_paths

    # ── Config ───────────────────────────────────────────────────────────────
    PATHS_VERSION  = "1"
    USE_BYPRODUCTS = True
    MAX_ATTEMPTS   = 10
    # ─────────────────────────────────────────────────────────────────────────

    all_paths = load_minpaths(f"paths_pv{PATHS_VERSION}")

    print("Building cross-feeding pair...")
    crossPair = build_crossfeeding_pair_from_paths(
        all_paths, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        use_byproducts=USE_BYPRODUCTS,
        max_attempts=MAX_ATTEMPTS)

    if crossPair is None:
        print("Failed to construct cross-feeding pair.")
    else:
        print(f"\nPair built successfully.")
        print(f"  cross_A: {len(crossPair['cross_A'])} rxns | cross_B: {len(crossPair['cross_B'])} rxns")
        print(f"  A donates met {crossPair['A_donated']} ({inv_met_map[crossPair['A_donated']]}) "
              f"--> B uses for core {crossPair['B_ext_core']} ({inv_met_map[crossPair['B_ext_core']]})")
        print(f"  B donates met {crossPair['B_donated']} ({inv_met_map[crossPair['B_donated']]}) "
              f"--> A uses for core {crossPair['A_ext_core']} ({inv_met_map[crossPair['A_ext_core']]})")

        t0 = time.time()
        result = splitByDemand_crossfeeding(
            stoich_matrix, rxnMat, prodMat, sumRxnVec,
            rho, pi, nutrientSet, Energy, Currency, Core, crossPair)
        t1 = time.time()

        print("\nYield results:")
        print(f"  Pair viable: {result['pair_viable']}")
        print(f"  Network A — viable: {result['viable_A']}, E: {result['E_A']:.6f}, B: {result['B_A']:.6f}")
        print(f"  Network B — viable: {result['viable_B']}, E: {result['E_B']:.6f}, B: {result['B_B']:.6f}")
        print(f"  Flux A→B: {result['flux_A_to_B']:.6f}, Flux B→A: {result['flux_B_to_A']:.6f}")
        print(f"  Time: {t1 - t0:.3f} s")
