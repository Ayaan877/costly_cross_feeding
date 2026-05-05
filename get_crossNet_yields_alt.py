from calculate_crossNet_yield_alt import *
from load_data import *
import sys
import time
import multiprocessing as mp
import pickle
import os
from directory_paths import resolve_crossnet_path, resolve_yield_path

'''
Calculates energy and biomass yields for a set of cross-feeding networks using the
original splitByDemand_alt method (greedy flux allocation, false viability check, 
non-limiting reactant check), for comparison with the edited version.
Optimized for parallel processing.
'''


def compute_crossfeeding_yield_alt(crossPair):
    result = splitByDemand_crossfeeding_alt(
        stoich_matrix, rxnMat, prodMat,
        sumRxnVec, rho, pi, nutrientSet,
        Energy, Currency, Core, crossPair)
    return (
        result['E_A'], result['B_A'], result['viable_A'],
        result['E_B'], result['B_B'], result['viable_B'],
        result['pair_viable'],
        result['flux_A_to_B'], result['flux_B_to_A'],
    )


if __name__ == "__main__":

    # Usage: get_alt_crossNet_yields.py <autonet_subdir> <autonet_file> <crossnet_subdir> <crossnet_file> <num_workers>
    #
    #   autonet_subdir  : "autonets_{source}_av{av}"
    #   autonet_file    : "{P|NP}_pv{pv}"  or  "P"
    #   crossnet_subdir : "crossnets_{source}_cv{cv}"
    #   crossnet_file   : "{byp|int}_{P|NP}"
    #   num_workers     : parallel worker count
    autonet_subdir  = sys.argv[1]
    autonet_file    = sys.argv[2]
    crossnet_subdir = sys.argv[3]
    crossnet_file   = sys.argv[4]
    num_workers     = int(sys.argv[5])

    crossnet_path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    output_path   = resolve_yield_path(autonet_subdir, autonet_file, "alt",
                                       crossnet_subdir, crossnet_file)
    os.makedirs(output_path.parent, exist_ok=True)

    with open(crossnet_path, "rb") as f:
        CrossNets = pickle.load(f)

    num_pairs = len(CrossNets)
    print(f"Loaded {num_pairs} cross-feeding pairs from {crossnet_path}")

    E_A_yields  = np.zeros(num_pairs)
    B_A_yields  = np.zeros(num_pairs)
    viable_A    = np.zeros(num_pairs, dtype=bool)
    E_B_yields  = np.zeros(num_pairs)
    B_B_yields  = np.zeros(num_pairs)
    viable_B    = np.zeros(num_pairs, dtype=bool)
    pair_viable = np.zeros(num_pairs, dtype=bool)
    flux_A_to_B = np.zeros(num_pairs)
    flux_B_to_A = np.zeros(num_pairs)

    start = time.time()

    print(f"Using {num_workers} parallel workers")
    with mp.Pool(processes=num_workers) as pool:
        for i, (EA, BA, vA, EB, BB, vB, vPair, fAB, fBA) in enumerate(
                pool.imap(compute_crossfeeding_yield_alt, CrossNets, chunksize=64)):
            E_A_yields[i]  = EA
            B_A_yields[i]  = BA
            viable_A[i]    = vA
            E_B_yields[i]  = EB
            B_B_yields[i]  = BB
            viable_B[i]    = vB
            pair_viable[i] = vPair
            flux_A_to_B[i] = fAB
            flux_B_to_A[i] = fBA

            if (i + 1) % 500 == 0:
                processed_ratio = (i + 1) / num_pairs
                viable_ratio    = np.sum(pair_viable[:i+1]) / (i + 1)
                print(f"  Processed {i + 1}/{num_pairs} ({processed_ratio:.2%}), "
                      f"pair viable: {np.sum(pair_viable[:i+1])}/{i + 1} ({viable_ratio:.2%})")

    elapsed = time.time() - start
    valid   = np.sum(pair_viable)
    print(f"\nCompleted in {elapsed:.2f} seconds")
    print(f"Viable pairs (both organisms produce all precursors): {valid}/{num_pairs}")

    with open(output_path, "wb") as f:
        pickle.dump({
            'E_A':         E_A_yields,
            'B_A':         B_A_yields,
            'viable_A':    viable_A,
            'E_B':         E_B_yields,
            'B_B':         B_B_yields,
            'viable_B':    viable_B,
            'pair_viable': pair_viable,
            'flux_A_to_B': flux_A_to_B,
            'flux_B_to_A': flux_B_to_A,
        }, f)
    print(f"Saved yields to {output_path}")
