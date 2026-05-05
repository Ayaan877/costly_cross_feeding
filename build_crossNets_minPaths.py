import sys
import time
from load_data import *
from datetime import datetime
from load_networks import load_minpaths
from generate_crossNets_minPaths import generate_crossNets_minPaths
from directory_paths import parse_crossnet_spec, resolve_crossnet_path

'''
Builder script for non-minimal cross-feeding networks from minimal core pathways.
Generates a specified number of unique cross-feeding pairs from a given set of
core minimal pathways, exchanging either byproduct or intermediate metabolites.
'''

if __name__ == "__main__":
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Usage: build_crossNets_minPaths.py <paths_subdir> <crossnet_subdir> <crossnet_file> <n_target> <n_workers>
    #
    #   paths_subdir    : "paths_pv{pv}"               e.g. "paths_pv2"
    #   crossnet_subdir : "crossnets_{source}_cv{cv}"  e.g. "crossnets_mp_cv1"
    #   crossnet_file   : "{byp|int}_NP"               (always NP: no post-assembly pruning)
    #   n_target        : number of unique cross-feeding pairs to generate
    #   n_workers       : number of parallel workers
    paths_subdir    = sys.argv[1]
    crossnet_subdir = sys.argv[2]
    crossnet_file   = sys.argv[3]
    n_target        = int(sys.argv[4])
    n_workers       = int(sys.argv[5])

    _, _, cross_type, pruning = parse_crossnet_spec(crossnet_subdir, crossnet_file)
    if cross_type not in ("byp", "int"):
        raise ValueError(f"Invalid cross_type in crossnet_file: '{cross_type}'. Must be 'byp' or 'int'.")
    if pruning != "NP":
        raise ValueError("build_crossNets_minPaths always generates NP (no post-assembly pruning) crossnets.")

    output_path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading MinPaths from {paths_subdir}...")
    all_paths = load_minpaths(paths_subdir)
    print(f"  Loaded {len(all_paths)} core targets, "
          f"{sum(len(p) for p in all_paths)} total pathways.")
    print(f"Output: {output_path}")

    start_time = time.time()

    pairs = generate_crossNets_minPaths(
        all_paths, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        n_target=n_target,
        n_workers=n_workers,
        save_path=output_path,
        use_byproducts=(cross_type == "byp"))

    total_time = time.time() - start_time
    print(f"Final: {len(pairs)} unique cross-feeding pairs generated")
    print(f"Total time: {total_time/60:.2f} minutes")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
