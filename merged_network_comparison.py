import pickle
import numpy as np
import os
from directory_paths import *

'''
Compare yields of cross-feeding pairs to their merged networks. Calculates the 
difference in biomass and energy yields between the sum of both cross-feeders
and the corresponding merged network. Identifies cases where the cross-feeding pair
outperforms the merged network, and quantifies how many pairs show this behavior. 
'''

AUTONET_SUBDIR  = "autonets_mp_av1"
AUTONET_FILE    = "NP_pv2"
CROSSNET_SUBDIR = "crossnets_mp_cv1"
CROSSNET_FILES  = ("byp_NP", "int_NP")   # both crossnet variants
YIELD_MODE      = "sbd"                # sbd | iter

results = {}

for crossnet_file in CROSSNET_FILES:
    _, cv, ctype, _ = parse_crossnet_spec(CROSSNET_SUBDIR, crossnet_file)

    crossnet_path     = resolve_crossnet_path(CROSSNET_SUBDIR, crossnet_file)
    cross_yield_path  = resolve_yield_path(AUTONET_SUBDIR, AUTONET_FILE, YIELD_MODE,
                                           CROSSNET_SUBDIR, crossnet_file)
    merged_yield_path = resolve_merged_yield_path(CROSSNET_SUBDIR, crossnet_file, YIELD_MODE)

    print(f"\n--- Processing CrossNet type '{ctype}' ---")
    print(f"  CrossNets        : {crossnet_path}")
    print(f"  Cross yields     : {cross_yield_path}")
    print(f"  Merged yields    : {merged_yield_path}")

    with open(crossnet_path, "rb") as f:
        CrossNets = pickle.load(f)

    with open(cross_yield_path, "rb") as f:
        cy = pickle.load(f)

    with open(merged_yield_path, "rb") as f:
        merged_E, merged_B, merged_viable = pickle.load(f)

    cross_E_A   = cy['E_A']
    cross_B_A   = cy['B_A']
    cross_E_B   = cy['E_B']
    cross_B_B   = cy['B_B']
    pair_viable = cy['pair_viable']

    n_pairs = len(CrossNets)
    print(f"  Pairs loaded    : {n_pairs}")
    print(f"  Pair-viable     : {int(np.sum(pair_viable))}")
    print(f"  Merged-viable   : {int(np.sum(merged_viable))}")

    # Valid: cross-feeding pair is viable AND merged network is viable
    valid_mask = pair_viable & merged_viable
    n_valid = int(np.sum(valid_mask))
    print(f"  Both viable     : {n_valid}")

    # Yield deltas: total pair yield minus merged network yield
    # Positive delta means the cross-feeding pair exceeds the merged network
    delta_B = np.full(n_pairs, np.nan)
    delta_E = np.full(n_pairs, np.nan)

    delta_B[valid_mask] = (cross_B_A[valid_mask] + cross_B_B[valid_mask]) - merged_B[valid_mask]
    delta_E[valid_mask] = (cross_E_A[valid_mask] + cross_E_B[valid_mask]) - merged_E[valid_mask]

    valid_indices = np.where(valid_mask)[0]
    vm = valid_mask

    outperformers_biomass = valid_indices[delta_B[vm] > 0]
    outperformers_energy  = valid_indices[delta_E[vm] > 0]
    outperformers_either  = valid_indices[(delta_B[vm] > 0) | (delta_E[vm] > 0)]

    print(f"  Pair outperforms merged (biomass): {len(outperformers_biomass)} / {n_valid} "
          f"({100*len(outperformers_biomass)/max(n_valid,1):.1f}%)")
    print(f"  Pair outperforms merged (energy) : {len(outperformers_energy)} / {n_valid} "
          f"({100*len(outperformers_energy)/max(n_valid,1):.1f}%)")
    print(f"  Pair outperforms merged (either) : {len(outperformers_either)} / {n_valid} "
          f"({100*len(outperformers_either)/max(n_valid,1):.1f}%)")

    results[ctype] = {
        'outperformers_biomass': outperformers_biomass,
        'outperformers_energy' : outperformers_energy,
        'outperformers_either' : outperformers_either,
        'delta_B'              : delta_B,
        'delta_E'              : delta_E,
        'merged_B'             : merged_B,
        'merged_E'             : merged_E,
        'pair_viable'          : pair_viable,
        'merged_viable'        : merged_viable,
        'valid_mask'           : valid_mask,
    }

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
source, av, pruning, _ = parse_autonet_spec(AUTONET_SUBDIR, AUTONET_FILE)
_, cv, _, _            = parse_crossnet_spec(CROSSNET_SUBDIR, CROSSNET_FILES[0])
save_dir  = YIELDS_DIR / f"yields_{source}_{YIELD_MODE}"
os.makedirs(save_dir, exist_ok=True)
save_path = save_dir / f"comparison_vs_merged_{pruning}_av{av}_cv{cv}.pkl"

with open(save_path, "wb") as f:
    pickle.dump(results, f)
print(f"\nSaved comparison data to {save_path}")
