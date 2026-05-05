from datetime import datetime
from reverse_scope import giveRevScope
from multiprocessing import Pool
import numpy as np
import pickle
import time

worker_data = {}

def init_worker(data):
    '''
    The large shared arrays (rxnMat, prodMat, etc.) and the pre-computed reverse scope
    (satRxns) are loaded once per worker process via init_worker to avoid
    redundant serialization. Only the RNG seed is passed per call, so each
    worker independently prunes the same reverse-scope subgraph to a different
    minimal network, yielding stochastic diversity across parallel calls.
    '''
    global worker_data
    worker_data = data


def single_variant(seed):
    '''
    Multiprocessing worker that generates one minimal pathway variant.
    '''
    d = worker_data
    rng = np.random.default_rng(seed)
    return d['randMinNetwork'](d['satRxns'], d['rxnMat'], d['prodMat'], d['sumRxnVec'],
                               d['target'], d['nutrientSet'], d['Currency'], rng=rng)


def generate_pruned_networks(target, rxnMat, prodMat, sumRxnVec,
                             nutrientSet, Currency, n_workers, randMinNetwork,
                             save_path=None, max_attempts=500,
                             plateau_window=5, plateau_threshold=2):
    '''
    Generates a set of unique minimal pathways from the nutrient set
    to a target core metabolite using parallel stochastic pruning. First
    computes the reverse scope of the target (giveRevScope), then dispatches
    n_workers parallel pruning calls per round with distinct seeds.
    Newly discovered networks are deduplicated by sorted reaction-tuple
    membership. Iteration stops when max_attempts rounds are completed or when
    the plateau detector finds fewer than plateau_threshold new unique networks
    in the last plateau_window rounds.

    target is a single metabolite index or array of indices. randMinNetwork is
    the pruning function to use (e.g. from batch_pruning or single_pruning).

    Returns a dict with keys: networks (list of numpy arrays), attempts (list
    of attempt counts), unique_counts (list of cumulative unique network counts).
    '''

    satMets, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec,

                                    nutrientSet, Currency, target)

    worker_data = dict(
        satRxns=satRxns, rxnMat=rxnMat, prodMat=prodMat,
        sumRxnVec=sumRxnVec, target=target,
        nutrientSet=nutrientSet, Currency=Currency,
        randMinNetwork=randMinNetwork,
    )

    unique_nets = set()

    attempts_list = []
    unique_counts = []

    attempt = 0

    with Pool(processes=n_workers, initializer=init_worker, initargs=(worker_data,)) as pool:
        while attempt < max_attempts:

            attempt += 1
            seeds = np.random.randint(0, 10**9, size=n_workers)

            start = time.time()

            new_nets = pool.map(single_variant, seeds)

            elapsed = time.time() - start

            for net in new_nets:
                unique_nets.add(tuple(sorted(net)))

            current_count = len(unique_nets)

            attempts_list.append(attempt)
            unique_counts.append(current_count)

            print(f"[{datetime.now().strftime('%H:%M:%S')}] "
                  f"Attempt {attempt}: {elapsed:.4f}s - "
                  f"{current_count} unique networks")

            if save_path is not None:
                results = {"networks": [np.array(net) for net in unique_nets],
                           "attempts": attempts_list,
                           "unique_counts": unique_counts}
                with open(save_path, "wb") as f:
                    pickle.dump(results, f)

            if len(unique_counts) > plateau_window:
                recent_growth = unique_counts[-1] - unique_counts[-plateau_window]
                if recent_growth <= plateau_threshold:
                    print("Unique network discovery has plateaued. Stopping...")
                    break

    return {"networks": [np.array(net) for net in unique_nets], 
            "attempts": attempts_list, 
            "unique_counts": unique_counts}