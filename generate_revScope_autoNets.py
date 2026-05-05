import numpy as np
import pickle
import time
from datetime import datetime
from multiprocessing import Pool
from reverse_scope import giveRevScope
from batch_pruning import randMinNetwork

worker_data = {}

def init_worker(data):
    global worker_data
    worker_data = data

def prune_worker(seed):
    '''
    Multiprocessing worker that generates an autonomous network by 
    pruning the reverse scope with a given seed.
    '''
    d = worker_data
    rng = np.random.default_rng(seed)
    return randMinNetwork(d['satRxns'], d['rxnMat'], d['prodMat'], d['sumRxnVec'],
                          d['coreTBPs'], d['nutrientSet'], d['Currency'], rng=rng)

def generate_revScopeAutoNets(rxnMat, prodMat, sumRxnVec, nutrientSet,
                              Currency, coreTBPs, n_target=50000,
                              n_workers=32, batch_size=None,
                              save_path=None, save_interval=1000):
    '''
    Generate `n_target` unique minimal autonomous networks by repeatedly
    batch-pruning the same reverse-scope subgraph with different seeds.
    '''
    if batch_size is None:
        batch_size = n_workers

    satMets, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec,
                                    nutrientSet, Currency, coreTBPs)
    print(f"Reverse scope found {int(np.sum(satRxns))} reactions.")

    unique_nets = set()
    attempts = 0
    processed = 0

    data = dict(satRxns=satRxns, rxnMat=rxnMat, prodMat=prodMat, sumRxnVec=sumRxnVec,
                coreTBPs=coreTBPs, nutrientSet=nutrientSet, Currency=Currency)

    with Pool(processes=n_workers, initializer=init_worker, initargs=(data,)) as pool:
        while len(unique_nets) < n_target:
            attempts += 1
            seeds = [int(s) for s in np.random.randint(0, 2**31, size=batch_size)]

            start = time.time()
            results = pool.map(prune_worker, seeds)
            elapsed = time.time() - start

            for net in results:
                unique_nets.add(tuple(sorted(net)))

            processed += len(results)
            print(f"[{datetime.now().strftime('%H:%M:%S')}] "
                  f"[Attempt {attempts}] {len(unique_nets)} unique autonets, {elapsed:.1f}s ")

            if save_path is not None and processed % save_interval < batch_size:
                with open(save_path, "wb") as f:
                    pickle.dump([np.array(net) for net in unique_nets], f)
                print(f"Checkpoint: {processed} processed, {len(unique_nets)} unique networks saved")

    networks = [np.array(net) for net in unique_nets]
    print(f"Finished: {len(networks)} unique autonets in {attempts} attempts.")

    if save_path is not None:
        with open(save_path, "wb") as f:
            pickle.dump(networks, f)
        print(f"Saved to {save_path}")
    return networks

