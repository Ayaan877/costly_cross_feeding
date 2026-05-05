from calculate_autoNet_yield import splitByDemand
from calculate_autoNet_yield_iterative import splitByDemandIterative
from load_data import *
import sys
import time
import multiprocessing as mp
import pickle
import os
from directory_paths import resolve_crossnet_path, resolve_merged_yield_path


def compute_yield_sbd(args):
    merged_net, nutrient_supply = args
    return splitByDemand(
        stoich_matrix, rxnMat, prodMat,
        sumRxnVec, rho, pi, nutrientSet,
        Energy, Currency, Core, merged_net,
        nutrientSupply=nutrient_supply)


def compute_yield_iter(args):
    merged_net, nutrient_supply = args
    return splitByDemandIterative(
        stoich_matrix, rxnMat, prodMat,
        sumRxnVec, rho, pi, nutrientSet,
        Energy, Currency, Core, merged_net,
        nutrientSupply=nutrient_supply)


def merge_pair(crossPair):
    """Return the union of both organisms' reaction indices and a nutrient
    supply dict that gives 2× supply to nutrients consumed by both organisms."""
    cross_A = np.array(crossPair['cross_A'], dtype=int)
    cross_B = np.array(crossPair['cross_B'], dtype=int)

    # Nutrients consumed as reactants by each organism's reactions.
    nuts_A = set(n for n in nutrientSet if rxnMat[cross_A][:, n].any())
    nuts_B = set(n for n in nutrientSet if rxnMat[cross_B][:, n].any())

    # Shared nutrients get double supply so total input matches the pair.
    nutrient_supply = {n: 2.0 for n in nuts_A & nuts_B}

    return np.union1d(cross_A, cross_B), nutrient_supply


if __name__ == "__main__":

    # Usage: get_merged_crossNet_yields.py <crossnet_subdir> <crossnet_file> <yield_mode> <num_workers>
    #
    #   crossnet_subdir : "crossnets_{source}_cv{cv}"   e.g. "crossnets_mp_cv1"
    #   crossnet_file   : "{byp|int}_{P|NP}"            e.g. "byp_P"
    #   yield_mode      : sbd | iter
    #   num_workers     : parallel worker count

    crossnet_subdir = "crossnets_mp_cv1" #sys.argv[1]
    crossnet_file   = "int_P" #sys.argv[2]
    yield_mode      = "sbd" #sys.argv[3]
    num_workers     = int(6) #sys.argv[4]

    if yield_mode == "sbd":
        compute_yield = compute_yield_sbd
    elif yield_mode == "iter":
        compute_yield = compute_yield_iter
    else:
        raise ValueError(f"Unknown yield_mode '{yield_mode}'. Use 'sbd' or 'iter'.")

    crossnet_path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    output_path   = resolve_merged_yield_path(crossnet_subdir, crossnet_file, yield_mode)
    os.makedirs(output_path.parent, exist_ok=True)

    with open(crossnet_path, "rb") as f:
        CrossNets = pickle.load(f)

    num_pairs = len(CrossNets)
    print(f"Loaded {num_pairs} cross-feeding pairs from {crossnet_path}")

    # Pre-compute all merged networks (union of both organisms' reactions).
    merged_net_args = [merge_pair(p) for p in CrossNets]

    E_yields  = np.zeros(num_pairs)
    B_yields  = np.zeros(num_pairs)
    viability = np.zeros(num_pairs, dtype=bool)

    start = time.time()

    print(f"Using {num_workers} parallel workers ({yield_mode} mode)")
    with mp.Pool(processes=num_workers) as pool:
        for i, (E_yield, B_yield, status) in enumerate(
                pool.imap(compute_yield, merged_net_args, chunksize=64)):
            E_yields[i]  = E_yield
            B_yields[i]  = B_yield
            viability[i] = status

            if (i + 1) % 500 == 0:
                processed_ratio = (i + 1) / num_pairs
                viable_ratio    = np.sum(viability[:i+1]) / (i + 1)
                print(f"  Processed {i + 1}/{num_pairs} ({processed_ratio:.2%}), "
                      f"viable: {np.sum(viability[:i+1])}/{i + 1} ({viable_ratio:.2%})")

    elapsed = time.time() - start
    valid = np.sum(viability)
    print(f"\nCompleted in {elapsed:.2f} seconds")
    print(f"Viable merged networks: {valid}/{num_pairs}")

    with open(output_path, "wb") as f:
        pickle.dump((E_yields, B_yields, viability), f)
    print(f"Saved merged yields to {output_path}")
