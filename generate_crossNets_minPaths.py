import numpy as np
import pickle
import time
import io
from contextlib import redirect_stdout
from datetime import datetime
from multiprocessing import Pool
from crossfeeding_minPaths import build_crossfeeding_pair_from_paths

worker_data = {}

def init_worker(data):
    global worker_data
    worker_data = data


def crossfeed_worker(seed):
    '''
    Worker function to generate a cross-feeding pair from random minimal core-producing pathways.
    '''
    d = worker_data
    np.random.seed(seed)

    with redirect_stdout(io.StringIO()):
        result = build_crossfeeding_pair_from_paths(
            d['all_paths'],
            d['rxnMat'], d['prodMat'], d['sumRxnVec'],
            d['nutrientSet'], d['Currency'], d['Core'],
            use_byproducts=d['use_byproducts'],
            max_attempts=d['max_attempts'])

    if result is None:
        return {'success': False, 'seed': seed}

    return {'success': True, 'seed': seed, 'result': result}


def pair_key(result):
    '''
    Canonical unordered key for deduplication by network content.
    '''
    return tuple(sorted([tuple(sorted(result['cross_A'])), tuple(sorted(result['cross_B']))]))


def generate_crossNets_minPaths(all_paths, rxnMat, prodMat, sumRxnVec,
                                nutrientSet, Currency, Core,
                                n_target, n_workers,
                                batch_size=None, save_path=None,
                                save_interval=1000, use_byproducts=False,
                                max_attempts=10):
    '''
    Generate up to `n_target` unique cross-feeding pairs from the library 
    of minimal core-producing pathways, using `n_workers` parallel processes.

    Each worker independently assembles two random autonomous networks by
    sampling one minimal pathway per core from the MinPaths library (8 pathways
    per organism) and attempts to construct a cross-feeding pair.
    '''
    if batch_size is None:
        batch_size = n_workers

    print(f"MinPaths library: {len(all_paths)} core targets, "
          f"{sum(len(p) for p in all_paths)} total pathways.")
    print(f"Target: {n_target} unique pairs | Workers: {n_workers} | "
          f"Batch: {batch_size} | use_byproducts: {use_byproducts}")

    unique_pairs = {}
    attempts = 0
    processed = 0
    total_successes = 0
    total_failures = 0

    rng = np.random.default_rng()

    data = dict(
        all_paths=all_paths,
        rxnMat=rxnMat, prodMat=prodMat, sumRxnVec=sumRxnVec,
        nutrientSet=nutrientSet, Currency=Currency, Core=Core,
        use_byproducts=use_byproducts,
        max_attempts=max_attempts,
    )

    with Pool(processes=n_workers, initializer=init_worker, initargs=(data,)) as pool:
        while len(unique_pairs) < n_target:
            attempts += 1

            seeds = rng.integers(0, 2**31, size=batch_size).tolist()

            start = time.time()
            results = pool.map(crossfeed_worker, seeds)
            elapsed = time.time() - start

            batch_new = 0
            for r in results:
                processed += 1
                if not r['success']:
                    total_failures += 1
                    print(f"  FAILED : attempt {processed}")
                else:
                    total_successes += 1
                    res = r['result']
                    key = pair_key(res)
                    is_new = key not in unique_pairs
                    if is_new:
                        unique_pairs[key] = res
                        batch_new += 1
                    print(f"  {'NEW   ' if is_new else 'DUP   '}: "
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
                _save(unique_pairs, save_path)
                print(f"  Checkpoint: {len(unique_pairs)} pairs saved to {save_path}")

    pairs = list(unique_pairs.values())
    print(f"\nFinished: {len(pairs)} unique pairs from {processed} attempts "
          f"({total_successes} successes, {total_failures} failures).")

    if save_path is not None:
        _save(unique_pairs, save_path)
        print(f"Final save: {len(pairs)} pairs -> {save_path}")

    return pairs


def _save(unique_pairs, save_path):
    pairs = list(unique_pairs.values())
    with open(save_path, "wb") as f:
        pickle.dump(pairs, f)
