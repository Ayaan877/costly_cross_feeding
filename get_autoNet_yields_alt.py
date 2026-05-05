from calculate_autoNet_yield_alt import *
from load_data import *
import sys
import time
import multiprocessing as mp
import pickle
import os
from directory_paths import resolve_autonet_path, resolve_yield_path

'''
Calculates energy and biomass yields for a set of autonomous networks using the
original splitByDemand_alt method (greedy flux allocation, false viability check, 
non-limiting reactant check), for comparison with the edited version. 
Optimized for parallel processing.
'''


def compute_yield_alt(net):
    return splitByDemand_alt(
        stoich_matrix, rxnMat, prodMat,
        sumRxnVec, rho, pi, nutrientSet,
        Energy, Currency, Core, net)


if __name__ == "__main__":

    # Usage: get_alt_autoNet_yields.py <autonet_subdir> <autonet_file> <num_workers>
    #
    #   autonet_subdir : "autonets_{source}_av{av}"
    #   autonet_file   : "{P|NP}_pv{pv}"  or  "P"
    #   num_workers    : parallel worker count
    autonet_subdir = sys.argv[1]
    autonet_file   = sys.argv[2]
    num_workers    = int(sys.argv[3])

    autonet_path = resolve_autonet_path(autonet_subdir, autonet_file)
    output_path  = resolve_yield_path(autonet_subdir, autonet_file, "alt")
    os.makedirs(output_path.parent, exist_ok=True)

    with open(autonet_path, "rb") as f:
        AutoNets = pickle.load(f)

    num_nets = len(AutoNets)
    print(f"Loaded {num_nets} networks from {autonet_path}")

    E_yields  = np.zeros(num_nets)
    B_yields  = np.zeros(num_nets)
    viability = np.zeros(num_nets, dtype=bool)

    start = time.time()

    print(f"Using {num_workers} parallel workers")
    with mp.Pool(processes=num_workers) as pool:
        for i, (E_yield, B_yield, status) in enumerate(
                pool.imap(compute_yield_alt, AutoNets, chunksize=64)):
            E_yields[i]  = E_yield
            B_yields[i]  = B_yield
            viability[i] = status

            if (i + 1) % 500 == 0:
                processed_ratio = (i + 1) / num_nets
                viable_ratio    = np.sum(viability[:i+1]) / (i + 1)
                print(f"  Processed {i + 1}/{num_nets} ({processed_ratio:.2%}), "
                      f"viable: {np.sum(viability[:i+1])}/{i + 1} ({viable_ratio:.2%})")

    elapsed = time.time() - start
    valid   = np.sum(viability)
    print(f"\nCompleted in {elapsed:.2f} seconds")
    print(f"Valid networks (all precursors produced): {valid}/{num_nets}")

    with open(output_path, "wb") as f:
        pickle.dump((E_yields, B_yields, viability), f)
    print(f"Saved yields to {output_path}")
