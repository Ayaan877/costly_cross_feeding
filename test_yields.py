import pickle
import time
from calculate_autoNet_yield import splitByDemand
from calculate_crossNet_yield import splitByDemand_crossfeeding
from directory_paths import resolve_autonet_path, resolve_crossnet_path
from load_data import *
from autonomy_check import verify_autonomy

'''
Tests yield calculations for autonomous and cross-feeding networks
'''

# ── Config ───────────────────────────────────────────────────────────────────
AUTONET_SUBDIR   = "autonets_rs_av1"   # autonets_{source}_av{version}
AUTONET_FILE     = "P"                  # P (rs is always pruned)
CROSSNET_SUBDIR  = "crossnets_rs_cv1"   # crossnets_{source}_cv{version}
CROSSNET_FILE    = "byp_P"              # {byp|int}_{P|NP}
AUTONET_IDX      = 49265
CROSSNET_IDX     = 1976
# ─────────────────────────────────────────────────────────────────────────────

autonet_path  = resolve_autonet_path(AUTONET_SUBDIR, AUTONET_FILE)
crossnet_path = resolve_crossnet_path(CROSSNET_SUBDIR, CROSSNET_FILE)

t1 = time.time()
with open(autonet_path, 'rb') as f:
    AutoNets = pickle.load(f)

autonet = AutoNets[AUTONET_IDX]

viable, missing_cores = verify_autonomy(autonet, rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)
produced_cores = [c for c in Core if c not in missing_cores]

print('\nSingle-network yield test')
print(f'Network index: {AUTONET_IDX}')

print(f'Autonomy check: {"All cores reachable" if viable else "All cores not reachable"}'
      f' | Missing cores: {missing_cores}')

E_single, B_single, status_single = splitByDemand(
    stoich_matrix, rxnMat, prodMat, sumRxnVec,
    rho, pi, nutrientSet, Energy, Currency, Core, autonet)

print(f'  Viable: {status_single}')
print(f'  E: {E_single}')
print(f'  B: {B_single}')
print(f'  Time taken: {time.time() - t1:.3f} s')

# Cross-feeding pair yield test
crossnet_idx = CROSSNET_IDX
t2 = time.time()
with open(crossnet_path, 'rb') as f:
    CrossNets = pickle.load(f)

crossPair = CrossNets[crossnet_idx]

viable_A, missing_A = verify_autonomy(crossPair['cross_A'], rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)
viable_B, missing_B = verify_autonomy(crossPair['cross_B'], rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, Core)

print('\nCross-feeding pair yield test')
print(f'Pair index: {crossnet_idx}')

print(f'Autonomy check A: {"All cores reachable" if viable_A else "All cores not reachable"} | Missing cores: {missing_A}')
print(f'Autonomy check B: {"All cores reachable" if viable_B else "All cores not reachable"} | Missing cores: {missing_B}')

result = splitByDemand_crossfeeding(
    stoich_matrix, rxnMat, prodMat, sumRxnVec,
    rho, pi, nutrientSet, Energy, Currency, Core, crossPair)

print(f'Network A size: {len(crossPair["cross_A"])} rxns, Network B size: {len(crossPair["cross_B"])} rxns')
print(f'Exchanged: A donates met {crossPair["A_donated"]} ({inv_met_map[crossPair["A_donated"]]}), '
      f'B donates met {crossPair["B_donated"]} ({inv_met_map[crossPair["B_donated"]]})')
print(f'  Pair viable: {result["pair_viable"]}')
print(f'  Network A - viable: {result["viable_A"]}, E: {result["E_A"]}, B: {result["B_A"]}')
print(f'  Network B - viable: {result["viable_B"]}, E: {result["E_B"]}, B: {result["B_B"]}')
print(f'  Flux A --> B: {result["flux_A_to_B"]:.6f}, Flux B --> A: {result["flux_B_to_A"]:.6f}')
print(f'  Time taken: {time.time() - t2:.3f} s')