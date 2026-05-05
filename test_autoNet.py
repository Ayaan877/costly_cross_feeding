'''
Builds a test minimal autonomous network from a set of minimal pathways
'''
from load_networks import load_minpaths
from combine_pathways import buildAutonomousNetwork
from load_data import *
import pickle
import time

if __name__ == "__main__":
    all_paths = load_minpaths("paths_pv2")
    paths = [all_paths[i][0] for i in range(8)]

    start = time.time()
    MinNet = buildAutonomousNetwork(paths, rxnMat, prodMat, sumRxnVec,
                                    nutrientSet, Currency, Core, prune=True)
    print(f"Time taken: {(time.time() - start)} seconds")

    with open('autoNet_test.pkl', "wb") as f:
            pickle.dump(MinNet, f)
