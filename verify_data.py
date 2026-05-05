from load_networks import *
from autonomy_check import verify_autonomy
from load_data import *

'''
Check Minimal Pathways
'''
PATHS_VERSION = "2"
file = "C00097"

with open(f"inv_rxn_map.pkl", "rb") as f:
    rxn_map = pickle.load(f)

with open(f"rxn_string_dict.pkl", "rb") as f:
    rxn_dict = pickle.load(f)

with open(resolve_paths_path(f"paths_pv{PATHS_VERSION}", file), "rb") as f:
    results = pickle.load(f)

min_paths = results['networks'] if isinstance(results, dict) else results

for path in min_paths[:5]:
    print(f"Pathway with {len(path)} reactions")
    for idx in path:
        rxn_id = rxn_map[idx]
        print(f"{rxn_id}: {rxn_dict[rxn_id]}")
    print('--------------------------------------')

'''
Check Crossfeeders
'''

dataset  = {"type": "cross", "subdir": "crossnets_rs_cv1", "file": "byp_P", "label": "Byproducts"}
PAIR_IDX = 17372

crossnet_path = resolve_crossnet_path(dataset['subdir'], dataset['file'])

with open(crossnet_path, "rb") as f:
    crossfeeders = pickle.load(f)

pair = crossfeeders[PAIR_IDX]

print(f"Cross-feeding pair index: {PAIR_IDX}")

# Check dependency
aug_A = nutrientSet + [pair['B_donated']]
aug_B = nutrientSet + [pair['A_donated']]

independent_A, dependent_A = verify_autonomy(pair['cross_A'], rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)
independent_B, dependent_B = verify_autonomy(pair['cross_B'], rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)

viability_A, missing_A = verify_autonomy(pair['cross_A'], rxnMat, prodMat, sumRxnVec, aug_A, Currency, Core)
viability_B, missing_B = verify_autonomy(pair['cross_B'], rxnMat, prodMat, sumRxnVec, aug_B, Currency, Core)

print(f"Pair viability: {bool(viability_A*viability_B)} | Missing cores: A={missing_A}, B={missing_B}")
print(f"Pair obligate dependence: {bool((not independent_A)*(not independent_B))} | Dependent cores: A={dependent_A}, B={dependent_B}")
print(f"    Network A size: auto {len(pair['auto_A'])} --> cross {len(pair['cross_A'])} reactions")
print(f"    Network B size: auto {len(pair['auto_B'])} --> cross {len(pair['cross_B'])} reactions")
print(f"    A needs metabolite {pair['B_donated']} from B to produce {pair['A_ext_core']}")
print(f"    B needs metabolite {pair['A_donated']} from A to produce {pair['B_ext_core']}")
print(f"    A --> B donor pathway has {len(pair['pathway_AB'])} reactions")
print(f"    B --> A donor pathway has {len(pair['pathway_BA'])} reactions")
print('--------------------------------------')

'''
Check Autonomous Networks
'''
dataset_auto = {"type": "auto", "subdir": "autonets_rs_av1", "file": "P", "label": "RS Minimal AutoNets"}
NET_IDX      = 17372

autonet_path = resolve_autonet_path(dataset_auto['subdir'], dataset_auto['file'])

with open(autonet_path, "rb") as f:
    autonets = pickle.load(f)

net = autonets[NET_IDX]

print(f"Autonomous network index: {NET_IDX}")
print(f"    Network size: {len(net)} reactions")

viable, missing = verify_autonomy(net, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)

print(f"    Viable: {viable} | Missing cores: {missing}")
print('--------------------------------------')