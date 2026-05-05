import sys
import time
from datetime import datetime
from generate_minPaths import generate_pruned_networks
from load_data import *
from batch_pruning import randMinNetwork
from directory_paths import resolve_paths_path

'''
Builder script for minimal core-producing pathways. Generates a specified number 
of unique minimal pathways for a given target metabolite.
'''

if __name__ == "__main__":
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Usage: build_minPaths.py <target_name> <paths_subdir> <n_workers> <max_attempts> <plateau_window> <plateau_threshold>
    #
    #   target_name  : KEGG metabolite ID  e.g. C00022
    #   paths_subdir : output subdirectory  e.g. paths_pv1
    target_name       = sys.argv[1]
    paths_subdir      = sys.argv[2]
    N_WORKERS         = int(sys.argv[3])
    max_attempts      = int(sys.argv[4])
    plateau_window    = int(sys.argv[5])
    plateau_threshold = int(sys.argv[6])

    target    = met_map[target_name]
    target_id = inv_met_map[target]

    output_path = resolve_paths_path(paths_subdir, target_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Running target: {target}")
    print(f"Target ID: {target_id}")
    print(f"Output: {output_path}")

    start_time = time.time()
    results = generate_pruned_networks(
        target, rxnMat, prodMat, sumRxnVec, nutrientSet,
        Currency, n_workers=N_WORKERS, randMinNetwork=randMinNetwork,
        save_path=output_path, max_attempts=max_attempts,
        plateau_window=plateau_window, plateau_threshold=plateau_threshold)

    if results:
        print(f"{len(results['networks'])} variants generated in {results['attempts'][-1]} attempts")

    total_time = time.time() - start_time
    print(f"Total time: {total_time/60:.2f} minutes")