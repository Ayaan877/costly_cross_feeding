import sys
import time
from datetime import datetime
from load_data import *
from directory_paths import parse_autonet_spec, resolve_autonet_path

'''
Builder script for autonomous networks. Generates a specified number of unique
autonomous networks using the revScope or minPath generation method, 
with optional post-assembly pruning.
'''

if __name__ == "__main__":
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Args supplied by PBS script (see run_autonomous_networks.pbs)
    # Usage: build_autoNets.py <autonet_subdir> <autonet_file> <n_target> <n_workers>
    #
    #   autonet_subdir : "autonets_{source}_av{version}"  e.g. "autonets_mp_av2"
    #   autonet_file   : "{P|NP}_pv{path_version}"  (mp)  or  "P"  (rs)
    #   n_target       : number of unique networks to generate
    #   n_workers      : parallel worker count
    autonet_subdir = sys.argv[1]
    autonet_file   = sys.argv[2]
    N_TARGET       = int(sys.argv[3])
    N_WORKERS      = int(sys.argv[4])

    source, av, pruning, pv = parse_autonet_spec(autonet_subdir, autonet_file)
    output_path = resolve_autonet_path(autonet_subdir, autonet_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    start_time = time.time()

    if source == "rs":
        from generate_revScope_autoNets import generate_revScopeAutoNets

        print(f"Generating revScope autonomous networks...")
        print(f"  Output: {output_path}")

        AutoNets = generate_revScopeAutoNets(
            rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core,
            n_target=N_TARGET, n_workers=N_WORKERS,
            save_path=output_path)

    else:  # mp
        from load_networks import load_minpaths
        from generate_minPath_autoNets import generate_minPathAutoNets

        paths_subdir = f"paths_pv{pv}"
        all_paths    = load_minpaths(paths_subdir)
        do_prune     = pruning == "P"

        print(f"Generating {'pruned' if do_prune else 'unpruned'} autonomous networks "
              f"from {paths_subdir}...")
        print(f"  Output: {output_path}")

        AutoNets = generate_minPathAutoNets(
            all_paths, rxnMat, prodMat, sumRxnVec,
            nutrientSet, Currency, Core, prune=do_prune,
            n_target=N_TARGET, n_workers=N_WORKERS,
            save_path=output_path)

    total_time = time.time() - start_time
    print(f"Saved {len(AutoNets)} autonomous networks to {output_path}")
    print(f"Total time: {total_time/60:.2f} minutes")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
