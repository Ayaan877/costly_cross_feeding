import numpy as np
import time
from prune_check import isCoreProduced

def randMinNetwork(satRxnVec, rxnMat, prodMat, sumRxnVec,
                   coreTBP, nutrientSet, Currency, rng=None):
    '''
    Prunes a set of satisfied reactions down to a minimal subgraph by
    repeatedly sweeping through all remaining reactions in a random order and
    removing each one that does not disrupt production of any target core
    metabolite. Sweeps continue until a full pass yields no removals, at which
    point the network is minimal under single-reaction removal.

    Returns a numpy array of reaction indices forming the minimal network.
    '''

    if rng is None:
        rng = np.random.default_rng()

    currSatRxnVec = np.copy(satRxnVec)
    print(f"Starting simple single pruning (Ayaan's version)...", flush=True)
    count = 0

    while True:
        count += 1
        # Keeping track of what the graph currently looks like.
        currSatRxns = np.nonzero(currSatRxnVec)[0]
        print(f"[Sweep {count}] Current size: {len(currSatRxns)} reactions", flush=True)

        removed_any = False

        # Randomize order for stochastic minimal networks
        removalOrder = rng.permutation(currSatRxns)

        n = len(removalOrder)
        print(f'[Sweep {count}] Checking reactions in random order for removability...', flush=True)
        removed = 0
        start = time.time()
        for i, remRxn in enumerate(removalOrder, 1):
            if isCoreProduced(remRxn, currSatRxnVec, rxnMat, prodMat, 
                              sumRxnVec, nutrientSet, Currency, coreTBP):

                currSatRxnVec[remRxn] = 0
                removed_any = True
                removed += 1

            if i == 1 or i % 100 == 0 or i == n:
                print(f"[Sweep {count}] {100*i/n:.1f}%, Time elapsed: {time.time() - start:.2f}s", flush=True)

        print(f"[Sweep {count}] Removed {removed} reactions this sweep.", flush=True)
        # If we completed a full pass with no removals → minimal
        if not removed_any:
            print("No more removable reactions → terminating.", flush=True)
            print(f'Minimal network size = {len(np.nonzero(currSatRxnVec)[0])}', flush=True)
            return np.nonzero(currSatRxnVec)[0]