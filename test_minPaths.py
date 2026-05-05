'''
Builds a test minimal core-producing network for a single target
'''
from reverse_scope import giveRevScope
from batch_pruning import randMinNetwork  
# Replace 'batch_pruning' with 'single_pruning' to test other methods
from load_data import *
import time

if __name__ == "__main__":

    with open("inv_rxn_map.pkl", "rb") as f:
        inv_rxn_map = pickle.load(f)

    target = Core[2] # Test on C00022 (Pyruvate)
    target_id = inv_met_map[target]
    target_name = cpd_string_dict[target_id]
    print(f"Testing pruning on target: {target_id} - {target_name}")
    
    start = time.time()
    satMets, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, target)
    print(f"Starting pruning...")
    print(f"Time taken for reverse scope: {time.time() - start:.2f} seconds")
    
    start = time.time()
    min_rxns = randMinNetwork(satRxns, rxnMat, prodMat, sumRxnVec, 
                             target, nutrientSet, Currency)
    print(f"Time taken for pruning: {time.time() - start:.2f} seconds")
    rxn_ids = [inv_rxn_map[idx] for idx in min_rxns]

    print("\nPrecursor Name:", target_id, target_name)
    print("Satisfied Metabolites:", int(np.sum(satMets)))
    print("Satisfied Reactions:", int(np.sum(satRxns)))
    print(f"Minimal Reactions ({len(rxn_ids)}): {rxn_ids}") 