import numpy as np
from satisfiability_check import markSatMetsRxns


def reachable_mets(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, excluded):
    '''
    Metabolites reachable from net via satisfiability check, minus excluded.
    '''
    net = np.asarray(net, dtype=int)
    if len(net) == 0:
        return set()
    satRxnVec = np.zeros(rxnMat.shape[0], dtype=int)
    satRxnVec[net] = 1
    satMets, _ = markSatMetsRxns(satRxnVec, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency)
    return set(np.nonzero(satMets)[0]) - set(excluded)


def get_byproducts(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, excluded,
                   context_rxns=None, other_rxns=None):
    '''
    Returns metabolites reachable from net (satisfiability check seeded from
    nutrientSet) that are not consumed by any reaction in context_rxns. If
    context_rxns is None, net itself is used as the consumption context, so
    the result is metabolites the network produces but never internally consumes.
    Metabolites in excluded are always filtered out.

    If other_rxns is provided, only metabolites also reachable from other_rxns
    are returned. This acts as a persistence check: the metabolite must have an
    independent source outside net and is therefore not exclusively produced by it.

    net and context_rxns may be a pathway (subset of reactions) or a full network.

    Returns a sorted numpy array of metabolite indices.
    '''
    net = np.asarray(net, dtype=int)
    produced = reachable_mets(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, excluded)

    ctx = net if context_rxns is None else np.asarray(context_rxns, dtype=int)
    consumed = (set(np.nonzero(np.logical_or.reduce(rxnMat[ctx]))[0])
                if len(ctx) > 0 else set())

    result = {m for m in produced if m not in consumed}

    if other_rxns is not None:
        other_rxns = np.asarray(other_rxns, dtype=int)
        persistent = reachable_mets(other_rxns, rxnMat, prodMat, sumRxnVec,
                                     nutrientSet, Currency, excluded)
        result = result & persistent

    return np.array(sorted(result), dtype=int)


def get_intermediates(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, excluded,
                      context_rxns=None, other_rxns=None):
    '''
    Returns metabolites reachable from net (satisfiability check seeded from
    nutrientSet) that are consumed by at least one reaction in context_rxns,
    thereby excluding pure byproducts. If context_rxns is None, net itself is
    used — returning metabolites that are both produced and internally consumed
    by the same reaction set. Metabolites in excluded are always filtered out.

    If other_rxns is provided, only metabolites also reachable from other_rxns
    are returned. This acts as a persistence check: the metabolite must have an
    independent source outside net and is therefore not exclusively produced by it.

    net and context_rxns may be a pathway (subset of reactions) or a full network.

    Returns a sorted numpy array of metabolite indices.
    '''
    net = np.asarray(net, dtype=int)
    produced = reachable_mets(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, excluded)

    ctx = net if context_rxns is None else np.asarray(context_rxns, dtype=int)
    consumed = (set(np.nonzero(np.logical_or.reduce(rxnMat[ctx]))[0])
                if len(ctx) > 0 else set())

    result = {m for m in produced if m in consumed}

    if other_rxns is not None:
        other_rxns = np.asarray(other_rxns, dtype=int)
        persistent = reachable_mets(other_rxns, rxnMat, prodMat, sumRxnVec,
                                     nutrientSet, Currency, excluded)
        result = result & persistent

    return np.array(sorted(result), dtype=int)


def get_candidates(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency,
                   Core, use_byproducts):
    '''
    Convenience wrapper returning candidate metabolites for cross-feeding
    secretion from a network. Currency, Core, and nutrientSet are always
    excluded. use_byproducts selects between strict byproducts (True: reachable
    but not internally consumed) and internal intermediates (False: reachable
    and internally consumed, so pure byproducts are excluded).

    Returns a numpy array of candidate metabolite indices.
    '''
    excluded = sorted(set(list(Currency) + list(Core) + list(nutrientSet)))

    if use_byproducts:
        return get_byproducts(net, rxnMat, prodMat, sumRxnVec,
                              nutrientSet, Currency, excluded)
    else:
        return get_intermediates(net, rxnMat, prodMat, sumRxnVec,
                                 nutrientSet, Currency, excluded)
    

# ----------------------------------------------
# Test Run
# ----------------------------------------------

if __name__ == "__main__":
    import pickle
    from load_data import *

    # ── Config ──────────────────────────────────────────────────────────────
    AUTONET_SUBDIR  = "autonets_rs_av1"   # autonets_{source}_av{version}
    AUTONET_FILE    = "P"                  # P (rs is always pruned)
    NET_IDX_A       = 2395
    NET_IDX_B       = 1965
    # ──────────────────────────────────────────────────────────────────
    from directory_paths import resolve_autonet_path
    with open(resolve_autonet_path(AUTONET_SUBDIR, AUTONET_FILE), "rb") as f:
        all_autonets = pickle.load(f)

    net_A = all_autonets[NET_IDX_A]
    net_B = all_autonets[NET_IDX_B]

    candidates_A = get_candidates(net_A, rxnMat, prodMat, sumRxnVec,
                                   nutrientSet, Currency, Core, use_byproducts=False)
    candidates_B = get_candidates(net_B, rxnMat, prodMat, sumRxnVec,
                                   nutrientSet, Currency, Core, use_byproducts=False)
    common = np.array(sorted(list(set(candidates_A) & set(candidates_B))), dtype=int)
    