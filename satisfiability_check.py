import numpy as np
from scipy import sparse

# Cache to avoid repeated dense→sparse conversion for the same arrays.
sparse_cache = {}

def make_sparse(arr):
    '''
    Convert a dense array to sparse format if not already sparse, with caching to
    optimize repeated calls on the same arrays.
    '''
    if sparse.issparse(arr):
        return arr
    key = id(arr)
    if key in sparse_cache and sparse_cache[key].shape == arr.shape:
        return sparse_cache[key]
    sp = sparse.csr_matrix(arr)
    sparse_cache[key] = sp
    return sp

def markSatMetsRxns(rxnProc, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency):
    '''
    Starting from a seed set of nutrients and currency metabolites, propagates
    reachability forward through the candidate reactions in rxnProc. A reaction
    becomes satisfied once ALL of its reactants are satisfied; a metabolite
    becomes satisfied once ANY reaction that produces it is satisfied. The
    procedure iterates until no new reactions or metabolites are marked.

    Accepts both dense (ndarray) and sparse (csr_matrix) inputs for rxnMat and
    prodMat; internally converts to sparse and compresses the metabolite space
    to only columns involved in active reactions, which speeds up large networks.

    rxnProc is a binary vector of length nRxns indicating the candidate reactions
    to consider (typically the full network or a pruned subset).

    Returns (satMetVec, satRxnVec), binary vectors of length nMets and nRxns
    respectively marking all reachable metabolites and reactions.
    '''
    n_rxns, n_mets = rxnMat.shape

    seeds = np.array(nutrientSet + Currency)

    # Only operate on candidate reactions to avoid full-matrix multiplies.
    active_rxns = np.nonzero(rxnProc)[0]

    if len(active_rxns) == 0:
        satMetVec = np.zeros(n_mets)
        satMetVec[seeds] = 1
        return satMetVec, np.zeros(n_rxns)

    # Convert to sparse once if needed (cached for repeated calls).
    sp_rxnMat = make_sparse(rxnMat)
    sp_prodMat = make_sparse(prodMat)

    # Extract submatrices for active reactions only.
    sub_rxnMat = sp_rxnMat[active_rxns]
    sub_prodMat = sp_prodMat[active_rxns]
    sub_sumRxnVec = sumRxnVec[active_rxns]

    # Compress metabolite dimension: keep only columns involved in any
    # active reaction plus the seed metabolites.
    involved_mets = np.union1d(
        np.union1d(sub_rxnMat.indices, sub_prodMat.indices), seeds)

    comp_rxnMat = sub_rxnMat[:, involved_mets]
    comp_prodMat_T = sub_prodMat[:, involved_mets].T.tocsr()
    
    # Map seed indices into compressed metabolite space.
    met_to_local = np.empty(n_mets, dtype=np.intp)
    met_to_local[involved_mets] = np.arange(len(involved_mets))

    n_involved = len(involved_mets)
    comp_satMetVec = np.zeros(n_involved)
    comp_satMetVec[met_to_local[seeds]] = 1

    comp_satRxnVec = np.zeros(len(active_rxns))

    while True:
        old_sub = comp_satRxnVec.copy()

        # Marking first reactions, then metabolites, iteratively.
        comp_satRxnVec = (np.asarray(comp_rxnMat.dot(comp_satMetVec)).ravel()
                          == sub_sumRxnVec) * 1.0
        comp_satMetVec = (np.asarray(comp_prodMat_T.dot(comp_satRxnVec)).ravel()
                          + comp_satMetVec > 0) * 1.0

        # Checking if all satisfied nodes have been marked.
        if np.array_equal(old_sub, comp_satRxnVec):
            break

    # Map back to full-size vectors.
    satMetVec = np.zeros(n_mets)
    satMetVec[involved_mets] = comp_satMetVec

    satRxnVec = np.zeros(n_rxns)
    satRxnVec[active_rxns] = comp_satRxnVec

    return satMetVec, satRxnVec