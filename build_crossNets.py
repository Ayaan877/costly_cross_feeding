import sys
import time
import pickle
from load_data import *
from datetime import datetime
from generate_crossNets import generate_crossNets
from directory_paths import parse_crossnet_spec, resolve_autonet_path, resolve_crossnet_path

'''
Builder script for minimal cross-feeding networks from autonomous networks. 
Generates a specified number of unique cross-feeding pairs from a given set of
autonomous networks, exchanging either byproduct or intermediate metabolites.
'''

if __name__ == "__main__":
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Usage: build_crossNets.py <autonet_subdir> <autonet_file> <crossnet_subdir> <crossnet_file> <n_target> <n_workers>
    #
    #   autonet_subdir  : "autonets_{source}_av{av}"  e.g. "autonets_mp_av1"
    #   autonet_file    : "{P|NP}_pv{pv}"  or  "P"
    #   crossnet_subdir : "crossnets_{source}_cv{cv}"  e.g. "crossnets_mp_cv1"
    #   crossnet_file   : "{byp|int}_{P|NP}"  e.g. "byp_NP"
    autonet_subdir  = sys.argv[1]
    autonet_file    = sys.argv[2]
    crossnet_subdir = sys.argv[3]
    crossnet_file   = sys.argv[4]
    n_target        = int(sys.argv[5])
    n_workers       = int(sys.argv[6])

    _, _, cross_type, _ = parse_crossnet_spec(crossnet_subdir, crossnet_file)
    if cross_type not in ("byp", "int"):
        raise ValueError(f"Invalid cross_type in crossnet_file: '{cross_type}'. Must be 'byp' or 'int'.")

    load_path   = resolve_autonet_path(autonet_subdir, autonet_file)
    output_path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(load_path, "rb") as f:
        all_autonets = pickle.load(f)

    print(f"Generating crossfeeder networks from {load_path}...")
    print(f"  Output: {output_path}")

    start_time = time.time()

    pairs = generate_crossNets(
        all_autonets, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        n_target=n_target,
        n_workers=n_workers,
        save_path=output_path,
        use_byproducts=(cross_type == "byp"))

    total_time = time.time() - start_time
    print(f"Final: {len(pairs)} unique cross-feeding pairs generated")
    print(f"Total time: {total_time/60:.2f} minutes")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

