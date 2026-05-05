"""
Determines the fraction of non-minimal cross-feeding pairs that are NOT obligate,
i.e. at least one organism in the pair is autonomous in the absence of the 
donated metabolite from its partner.

Uses verify_autonomy (reachability-based) consistent with how obligate
dependence is enforced during cross-feeding construction.
"""

import pickle
import numpy as np
from load_data import *
from autonomy_check import verify_autonomy
from directory_paths import resolve_crossnet_path

CROSSNET_SUBDIR = "crossnets_mp_cv1"
CROSSNET_FILES  = ("byp_NP", "int_NP")

for crossnet_file in CROSSNET_FILES:
    path = resolve_crossnet_path(CROSSNET_SUBDIR, crossnet_file)
    with open(path, "rb") as f:
        CrossNets = pickle.load(f)

    n_pairs     = len(CrossNets)
    n_A_auto    = 0   # A is autonomous without B_donated
    n_B_auto    = 0   # B is autonomous without A_donated
    n_not_oblig = 0   # at least one is autonomous (pair is not obligate)

    for pair in CrossNets:
        cross_A   = pair['cross_A']
        cross_B   = pair['cross_B']
        # A receives B_donated; check if A is viable WITHOUT it
        a_auto, _ = verify_autonomy(cross_A, rxnMat, prodMat, sumRxnVec,
                                    nutrientSet, Currency, Core)
        # B receives A_donated; check if B is viable WITHOUT it
        b_auto, _ = verify_autonomy(cross_B, rxnMat, prodMat, sumRxnVec,
                                    nutrientSet, Currency, Core)

        if a_auto:
            n_A_auto += 1
        if b_auto:
            n_B_auto += 1
        if a_auto or b_auto:
            n_not_oblig += 1

    frac_not_oblig = n_not_oblig / n_pairs if n_pairs > 0 else float('nan')
    frac_A_auto    = n_A_auto    / n_pairs if n_pairs > 0 else float('nan')
    frac_B_auto    = n_B_auto    / n_pairs if n_pairs > 0 else float('nan')

    print(f"\n=== {crossnet_file} ===")
    print(f"  Total pairs             : {n_pairs}")
    print(f"  A autonomous (no B_don) : {n_A_auto:>6}  ({100*frac_A_auto:.2f}%)")
    print(f"  B autonomous (no A_don) : {n_B_auto:>6}  ({100*frac_B_auto:.2f}%)")
    print(f"  Not obligate (either)   : {n_not_oblig:>6}  ({100*frac_not_oblig:.2f}%)")
    print(f"  Obligate (both depend)  : {n_pairs - n_not_oblig:>6}  "
          f"({100*(1 - frac_not_oblig):.2f}%)")
