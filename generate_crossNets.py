import numpy as np
import pickle
import time
import io
from contextlib import redirect_stdout
from datetime import datetime
from multiprocessing import Pool
from crossfeeding import build_crossfeeding_pair

worker_data = {}

def init_worker(data):
    global worker_data
    worker_data = data

def crossfeed_worker(args):
    '''
    Worker function to generate a cross-feeding pair from random minimal autonomous networks.
    '''
    i_A, i_B, seed = args
    d = worker_data
    np.random.seed(seed)

    net_A = d['autonets'][i_A]
    net_B = d['autonets'][i_B]

    with redirect_stdout(io.StringIO()):
        result = build_crossfeeding_pair(
            net_A, net_B,
            d['rxnMat'], d['prodMat'], d['sumRxnVec'],
            d['nutrientSet'], d['Currency'], d['Core'],
            use_byproducts=d['use_byproducts'], max_attempts=10, max_runs=3)

    if result is None:
        return {
            'success': False,
            'i_A': i_A,
            'i_B': i_B,
            'auto_A_size': len(net_A),
            'auto_B_size': len(net_B),}

    return {
        'success': True,
        'i_A': i_A,
        'i_B': i_B,
        'result': result,
    }


def pair_key(result):
    '''
    Canonical unordered key for deduplication by network content.
    '''
    return tuple(sorted([tuple(sorted(result['cross_A'])), tuple(sorted(result['cross_B']))]))


def generate_crossNets(autonets, rxnMat, prodMat, sumRxnVec,
                                nutrientSet, Currency, Core,
                                n_target, n_workers,
                                batch_size=None, save_path=None,
                                save_interval=1000, use_byproducts=False):
    '''
    Generate up to `n_target` unique cross-feeding pairs from the library 
    of minimal autonomous networks, using `n_workers` parallel processes.

    Each worker independently samples two random autonomous networks 
    and attempts to construct a minimal cross-feeding pair.
    '''
    if batch_size is None:
        batch_size = n_workers

    n_nets = len(autonets)
    print(f"Loaded {n_nets} autonomous networks.")
    print(f"Target: {n_target} unique pairs | Workers: {n_workers} | Batch: {batch_size}")

    unique_pairs = {}   # key -> result dict
    attempts = 0
    processed = 0
    total_successes = 0
    total_failures = 0

    rng = np.random.default_rng()

    data = dict(autonets=autonets, rxnMat=rxnMat, prodMat=prodMat,
                sumRxnVec=sumRxnVec, nutrientSet=nutrientSet,
                Currency=Currency, Core=Core, use_byproducts=use_byproducts)

    with Pool(processes=n_workers, initializer=init_worker, initargs=(data,)) as pool:
        while len(unique_pairs) < n_target:
            attempts += 1

            # Sample batch_size distinct ordered pairs (i_A != i_B)
            idx = rng.integers(0, n_nets, size=(batch_size, 2))
            # Re-roll any i_A == i_B
            same = idx[:, 0] == idx[:, 1]
            while same.any():
                idx[same] = rng.integers(0, n_nets, size=(same.sum(), 2))
                same = idx[:, 0] == idx[:, 1]

            seeds = rng.integers(0, 2**31, size=batch_size)
            job_args = [(int(idx[k, 0]), int(idx[k, 1]), int(seeds[k]))
                        for k in range(batch_size)]

            start = time.time()
            results = pool.map(crossfeed_worker, job_args)
            elapsed = time.time() - start

            batch_new = 0
            for r in results:
                processed += 1
                if not r['success']:
                    total_failures += 1
                    print(f"  FAILED : nets [{r['i_A']}, {r['i_B']}] "
                          f"(sizes {r['auto_A_size']}, {r['auto_B_size']})")
                else:
                    total_successes += 1
                    res = r['result']
                    key = pair_key(res)
                    is_new = key not in unique_pairs
                    if is_new:
                        unique_pairs[key] = res
                        batch_new += 1
                    print(f"  {'NEW   ' if is_new else 'DUP   '}: "
                          f"nets [{r['i_A']}, {r['i_B']}] | "
                          f"A: {len(res['auto_A'])} --> {len(res['cross_A'])} rxns | "
                          f"B: {len(res['auto_B'])} --> {len(res['cross_B'])} rxns | "
                          f"A donates met {res['A_donated']} | "
                          f"B donates met {res['B_donated']}")

            print(f"[{datetime.now().strftime('%H:%M:%S')}] "
                  f"Batch {attempts}: {batch_new} new | "
                  f"{len(unique_pairs)} unique / {processed} processed "
                  f"({total_successes} success, {total_failures} fail) | "
                  f"{elapsed:.1f}s")

            if (save_path is not None
                    and processed % save_interval < batch_size
                    and len(unique_pairs) > 0):
                save_data(unique_pairs, save_path)
                print(f"  Checkpoint: {len(unique_pairs)} pairs saved to {save_path}")

    pairs = list(unique_pairs.values())
    print(f"\nFinished: {len(pairs)} unique pairs from {processed} attempts "
          f"({total_successes} successes, {total_failures} failures).")

    if save_path is not None:
        save_data(unique_pairs, save_path)
        print(f"Final save: {len(pairs)} pairs -> {save_path}")

    return pairs


def save_data(unique_pairs, save_path):
    pairs = list(unique_pairs.values())
    with open(save_path, "wb") as f:
        pickle.dump(pairs, f)
