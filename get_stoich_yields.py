"""
Computes stoichiometric biomass (biomassCost) and energy (fitCost) yields for
autonomous networks and cross-feeding pairs using the raw stoichiometric
approach from fit_cost.py.

Biomass yield (biomassCost)
---------------------------
    numerator   = int( sum_{r} sum_{m in Core}  S[r, m] )
    denominator = |sum_{r: net nutrient uptake < 0}  sum_{m in nutrientSet} S[r, m]|
    yield       = numerator / denominator   (0.0 if denominator == 0)

Energy yield (fitCost)
----------------------
    numerator   = int( sum_{r} sum_{m in Energy}  S[r, m] * mulFac[m] )
    mulFac      = [1, 3, ...]  (ATP weight 1; redox carriers weight 3 per NADPH≈3ATP)
    denominator = same as above
    yield       = numerator / denominator

For cross-feeding pairs the donated metabolite is added to each organism's
effective nutrient set before computing both yields.

Outputs
-------
data/yields/stoich_auto_{autoID}.pkl
    dict: {biomass: ndarray(n_nets,), fitness: ndarray(n_nets,)}

data/yields/stoich_cross_{autoID}_{ctype}_v{CROSSNET_ID}.pkl
    dict: {biomass_A: ndarray, fitness_A: ndarray,
           biomass_B: ndarray, fitness_B: ndarray}
"""

import numpy as np
import pickle
import multiprocessing as mp
import time
import sys
import os

from load_data import stoich_matrix, nutrientSet, Energy, Core

# ---------------------------------------------------------------------------
# Stoichiometric yield functions (adapted from fit_cost.py)
# ---------------------------------------------------------------------------

# ATP contributes 1 energy-equivalent; redox carriers contribute 3 (NADPH ≈ 3 ATP).
MULTI_FAC = np.array([1, 3, 3], dtype=float)


def nutrients_consumed(rxns, nut_set):
    """
    Total absolute flux into nutrient metabolites across the network.
    Only reactions whose net stoichiometry on nut_set is negative (i.e. they
    consume nutrients) are counted.
    """
    return abs(sum(
        float(np.sum(stoich_matrix[rxn][nut_set]))
        for rxn in rxns
        if float(np.sum(stoich_matrix[rxn][nut_set])) < 0.0
    ))


def biomass_cost(rxns, nut_set=None):
    """
    Net core-metabolite production per unit nutrient consumed.

    Parameters
    ----------
    rxns    : array-like of int, reaction indices
    nut_set : list of int, metabolite indices to use as nutrients
              (defaults to global nutrientSet from load_data)
    """
    if nut_set is None:
        nut_set = nutrientSet
    rxns = np.asarray(rxns, dtype=int)
    denom = nutrients_consumed(rxns, nut_set)
    if denom == 0.0:
        return 0.0
    numer = int(sum(float(np.sum(stoich_matrix[rxn][Core])) for rxn in rxns))
    return numer / denom


def fit_cost(rxns, nut_set=None):
    """
    Net energy-equivalent production per unit nutrient consumed.

    Parameters
    ----------
    rxns    : array-like of int, reaction indices
    nut_set : list of int, metabolite indices to use as nutrients
    """
    if nut_set is None:
        nut_set = nutrientSet
    rxns = np.asarray(rxns, dtype=int)
    denom = nutrients_consumed(rxns, nut_set)
    if denom == 0.0:
        return 0.0
    numer = int(sum(
        float(np.dot(stoich_matrix[rxn][Energy], MULTI_FAC))
        for rxn in rxns
    ))
    return numer / denom


def comp_fitness(rxns, nut_set=None, alpha=1.0, beta=0.0, gamma=0.10):
    """
    Composite fitness score.

    alpha * fit_cost + beta * biomass_cost - gamma * |network|

    Default weights: alpha=1, beta=0, gamma=0.1.
    Network-size penalty (gamma) discourages unnecessarily large networks.
    """
    if nut_set is None:
        nut_set = nutrientSet
    rxns = np.asarray(rxns, dtype=int)
    return (alpha * fit_cost(rxns, nut_set)
            + beta * biomass_cost(rxns, nut_set)
            - gamma * len(rxns))


# ---------------------------------------------------------------------------
# Per-network worker functions — must be at module level for Windows pickling
# ---------------------------------------------------------------------------

def compute_auto_yields(net):
    """Return (biomass_cost, fit_cost) for a single autonomous network."""
    rxns = np.asarray(net, dtype=int)
    return biomass_cost(rxns), fit_cost(rxns)


def compute_cross_yields(crossPair):
    """
    Return (bms_A, fit_A, bms_B, fit_B) for a cross-feeding pair.

    Each organism's effective nutrient set is extended by the metabolite it
    receives from its partner.
    """
    cross_A   = np.asarray(crossPair['cross_A'], dtype=int)
    cross_B   = np.asarray(crossPair['cross_B'], dtype=int)
    B_donated = int(crossPair['B_donated'])   # metabolite A receives from B
    A_donated = int(crossPair['A_donated'])   # metabolite B receives from A

    nut_A = list(set(nutrientSet) | {B_donated})
    nut_B = list(set(nutrientSet) | {A_donated})

    return (
        biomass_cost(cross_A, nut_A),
        fit_cost(cross_A, nut_A),
        biomass_cost(cross_B, nut_B),
        fit_cost(cross_B, nut_B),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
# Usage (auto, rs):  get_stoich_yields.py auto rs  <autonet_id> <num_workers>
# Usage (auto, mp):  get_stoich_yields.py auto mp  <autonet_id> <num_workers> <pruner> <pruning> <paths_version>
# Usage (cross, rs): get_stoich_yields.py cross rs <autonet_id> <crossnet_id> <crossnet_type> <num_workers>
# Usage (cross, mp): get_stoich_yields.py cross mp <autonet_id> <crossnet_id> <crossnet_type> <num_workers> <pruner> <pruning> <paths_version>

if __name__ == "__main__":

    # Usage: get_stoich_yields.py <net_type> <autonet_subdir> <autonet_file> <num_workers> [<crossnet_subdir> <crossnet_file>]
    #
    #   net_type        : auto | cross
    #   autonet_subdir  : "autonets_{source}_av{av}"
    #   autonet_file    : "{P|NP}_pv{pv}"  or  "P"
    #   num_workers     : parallel worker count
    #   crossnet_subdir : "crossnets_{source}_cv{cv}"  (cross only)
    #   crossnet_file   : "{byp|int}_{P|NP}"            (cross only)
    from directory_paths import (resolve_autonet_path, resolve_crossnet_path,
                                 resolve_yield_path)

    net_type       = sys.argv[1]
    autonet_subdir = sys.argv[2]
    autonet_file   = sys.argv[3]
    NUM_WORKERS    = int(sys.argv[4])

    if net_type not in ("auto", "cross"):
        raise ValueError(f"Unknown net_type '{net_type}'. Use 'auto' or 'cross'.")

    if net_type == "cross":
        crossnet_subdir = sys.argv[5]
        crossnet_file   = sys.argv[6]

    autonet_path = resolve_autonet_path(autonet_subdir, autonet_file)
    auto_out     = resolve_yield_path(autonet_subdir, autonet_file, "stoich")
    os.makedirs(auto_out.parent, exist_ok=True)

    # -----------------------------------------------------------------------
    # AutoNets — stoichiometric yields + autonomous ceiling
    # -----------------------------------------------------------------------
    with open(autonet_path, "rb") as f:
        AutoNets = pickle.load(f)
    n_auto = len(AutoNets)
    print(f"Loaded {n_auto} networks from {autonet_path}")

    auto_biomass = np.zeros(n_auto)
    auto_fitness = np.zeros(n_auto)

    t0 = time.time()
    print(f"Using {NUM_WORKERS} parallel workers")
    with mp.Pool(processes=NUM_WORKERS) as pool:
        for i, (bms, fit) in enumerate(
                pool.imap(compute_auto_yields, AutoNets, chunksize=128)):
            auto_biomass[i] = bms
            auto_fitness[i] = fit
            if (i + 1) % 500 == 0:
                processed_ratio = (i + 1) / n_auto
                print(f"  Processed {i + 1}/{n_auto} ({processed_ratio:.2%})")

    elapsed = time.time() - t0
    max_biomass = float(np.max(auto_biomass))
    max_fitness = float(np.max(auto_fitness))
    print(f"\nCompleted in {elapsed:.2f} seconds")
    print(f"Autonomous ceiling — max biomass: {max_biomass:.6f}  |  max fitness: {max_fitness:.6f}")

    with open(auto_out, "wb") as f:
        pickle.dump({"biomass": auto_biomass, "fitness": auto_fitness}, f)
    print(f"Saved yields to {auto_out}")

    if net_type == "cross":
        # -------------------------------------------------------------------
        # CrossNets — compare each organism against the autonomous ceiling
        # -------------------------------------------------------------------
        cross_path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
        cross_out  = resolve_yield_path(autonet_subdir, autonet_file, "stoich",
                                        crossnet_subdir, crossnet_file)
        os.makedirs(cross_out.parent, exist_ok=True)

        with open(cross_path, "rb") as f:
            CrossNets = pickle.load(f)
        n_pairs = len(CrossNets)
        print(f"Loaded {n_pairs} cross-feeding pairs from {cross_path}")

        bms_A = np.zeros(n_pairs)
        fit_A = np.zeros(n_pairs)
        bms_B = np.zeros(n_pairs)
        fit_B = np.zeros(n_pairs)

        t0 = time.time()
        print(f"Using {NUM_WORKERS} parallel workers")
        with mp.Pool(processes=NUM_WORKERS) as pool:
            for i, (ba, fa, bb, fb) in enumerate(
                    pool.imap(compute_cross_yields, CrossNets, chunksize=64)):
                bms_A[i] = ba
                fit_A[i] = fa
                bms_B[i] = bb
                fit_B[i] = fb
                if (i + 1) % 500 == 0:
                    processed_ratio = (i + 1) / n_pairs
                    both_beat_so_far = int(np.sum((bms_A[:i+1] > max_biomass) & (bms_B[:i+1] > max_biomass)))
                    print(f"  Processed {i + 1}/{n_pairs} ({processed_ratio:.2%}), "
                          f"both beat ceiling: {both_beat_so_far}/{i + 1}")

        elapsed = time.time() - t0

        both_beat   = (bms_A > max_biomass) & (bms_B > max_biomass)
        either_beat = (bms_A > max_biomass) | (bms_B > max_biomass)
        print(f"\nCompleted in {elapsed:.2f} seconds")
        print(f"Both organisms beat ceiling  : {int(np.sum(both_beat))}/{n_pairs}")
        print(f"Either organism beats ceiling: {int(np.sum(either_beat))}/{n_pairs}")

        with open(cross_out, "wb") as f:
            pickle.dump({
                "biomass_A": bms_A, "fitness_A": fit_A,
                "biomass_B": bms_B, "fitness_B": fit_B,
            }, f)
        print(f"Saved yields to {cross_out}")
