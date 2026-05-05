import numpy as np
from reverse_scope import giveRevScope
from batch_pruning import randMinNetwork
from autonomy_check import verify_autonomy
from find_intermediates import get_byproducts


def make_donor_pathway(donor_met, target_core, rxnMat, prodMat, sumRxnVec, Currency):
    """
    Build a minimal pathway from donor_met to target_core using only
    `donor_met` (plus currency) as the nutrient source, via reverse scope
    expansion and random minimal-network pruning.
    """
    try:
        _, satRxns = giveRevScope(rxnMat, prodMat, sumRxnVec,
                                  [donor_met], Currency, target_core)
    except ValueError:
        return None

    if np.sum(satRxns) == 0:
        return None

    pathway = randMinNetwork(satRxns, rxnMat, prodMat, sumRxnVec,
                             target_core, [donor_met], Currency)

    # Verify the pathway actually consumes the donor metabolite.
    if not any(donor_met in set(np.nonzero(rxnMat[r])[0]) for r in pathway):
        print(f"  ERROR: Donor pathway does not consume met {donor_met}.", flush=True)
        return None

    return pathway


def randAutoNet(all_paths, Core):
    """
    Assemble a random autonomous network by one random pathway per core 
    and taking their union.
    """
    path_dict = {}
    for i, core in enumerate(Core):
        paths = all_paths[i]
        path_dict[core] = paths[np.random.randint(len(paths))]
    full_rxns = np.unique(np.concatenate(list(path_dict.values()))).astype(int)
    return path_dict, full_rxns


def union_paths(path_dict, Core, exclude_core=None):
    """Union of all pathway reaction sets, optionally excluding one core's pathway."""
    arrays = [path_dict[c] for c in Core if c != exclude_core]
    if not arrays:
        return np.array([], dtype=int)
    return np.unique(np.concatenate(arrays)).astype(int)


# ── Main construction function ─────────────────────────────────────────────────

def build_crossfeeding_pair_from_paths(all_paths, rxnMat, prodMat, sumRxnVec,
                                       nutrientSet, Currency, Core,
                                       use_byproducts=True,
                                       max_attempts=10):
    """
    Build crossfeeding pairs by replacing one core pathway in each organism 
    with a donor pathway (byproduct/intermediate --> core):

      Step 1: Build random net_A and net_B from minPaths.
      Step 2: Find int_A - a byproduct/intermediate from net_A to some core c_B.
      Step 3: Generate donor pathway int_A --> c_B and 
              replace c_B's pathway in net_B with it.
      Step 4: Find int_B - a byproduct/intermediate from donor_pathway_AB to some core c_A.
      Step 5: Generate donor pathway int_B --> c_A (c_A ≠ c_B) and 
              replace c_A's pathway in net_A with it.
      Step 6: Assemble:
                cross_A = A's 7 MinPaths (excl. c_A) + donor_path_BA
                cross_B = B's 7 MinPaths (excl. c_B) + donor_pathway_AB

    When use_byproducts=False (intermediate mode), the exchange metabolite must
    be produced by produced by at least one other pathway besides the core producing
    pathway being replaced.
    """
    Core = list(Core)
    excluded = sorted(set(list(Core) + list(nutrientSet) + list(Currency)))

    for attempt in range(1, max_attempts + 1):
        print(f"\n[Attempt {attempt}/{max_attempts}]", flush=True)

        # --- Step 1: Build two random organisms from MinPaths ---
        path_dict_A, full_rxns_A = randAutoNet(all_paths, Core)
        path_dict_B, full_rxns_B = randAutoNet(all_paths, Core)

        print(f"  net_A: {len(full_rxns_A)} rxns | net_B: {len(full_rxns_B)} rxns",
              flush=True)

        # --- Step 2: build donor_pathway_AB (int_A --> c_B) ---
        donor_pathway_AB = None
        int_A = None
        c_B = None

        for core_B in np.random.permutation(Core):
            int_A_origin_path = path_dict_A[core_B]

            if use_byproducts:
                int_A_candidates = get_byproducts(
                    int_A_origin_path, rxnMat, prodMat, sumRxnVec,
                    nutrientSet, Currency, excluded,
                    context_rxns=full_rxns_A)
            else:
                other_rxns_A = union_paths(path_dict_A, Core, exclude_core=core_B)
                int_A_candidates = get_byproducts(
                    int_A_origin_path, rxnMat, prodMat, sumRxnVec,
                    nutrientSet, Currency, excluded,
                    context_rxns=full_rxns_A, other_rxns=other_rxns_A)

            if len(int_A_candidates) == 0:
                continue

            for met in np.random.permutation(int_A_candidates):
                donor_path = make_donor_pathway(
                                int(met), core_B, rxnMat, prodMat, sumRxnVec, Currency)
                if donor_path is not None:
                    donor_pathway_AB = donor_path
                    int_A = int(met)
                    c_B = core_B
                    break

            if donor_pathway_AB is not None:
                break

        if donor_pathway_AB is None or c_B is None:
            print("  No donor pathway found for B.", flush=True)
            continue

        print(f"  Donor path B: int_A={int_A} --> c_B={c_B} ({len(donor_pathway_AB)} rxns)",
              flush=True)

        # Build modified net_B: 7 original paths (excl. c_B) + donor_pathway_AB.
        B_remaining_rxns = union_paths(path_dict_B, Core, exclude_core=c_B)
        cross_B = np.unique(np.concatenate([B_remaining_rxns, donor_pathway_AB])).astype(int)

        # --- Step 3: Find int_B from donor_pathway_AB and build donor_path_BA ---
        # Note: int_B is pulled from pathway_AB, not the entire cross_B
        if use_byproducts:
            int_B_candidates = get_byproducts(
                donor_pathway_AB, rxnMat, prodMat, sumRxnVec,
                nutrientSet, Currency, excluded,
                context_rxns=cross_B)
        else:
            # Persistence constraint: int_B must also be reachable from B's remaining
            # pathways (not exclusively from donor_pathway_AB).
            int_B_candidates = get_byproducts(
                donor_pathway_AB, rxnMat, prodMat, sumRxnVec,
                nutrientSet, Currency, excluded,
                context_rxns=cross_B, other_rxns=B_remaining_rxns)

        if len(int_B_candidates) == 0:
            print("  No int_B candidates from donor_pathway_AB.", flush=True)
            continue

        donor_path_BA = None
        int_B = None
        c_A = None

        for core_A in np.random.permutation([c for c in Core if c != c_B]):
            for met in np.random.permutation(int_B_candidates):
                donor_path = make_donor_pathway(
                                int(met), core_A, rxnMat, prodMat, sumRxnVec, Currency)
                if donor_path is not None:
                    donor_path_BA = donor_path
                    int_B = int(met)
                    c_A = core_A
                    break
            if donor_path_BA is not None:
                break

        if donor_path_BA is None or c_A is None:
            print("  No donor pathway found for A.", flush=True)
            continue

        print(f"  Donor path A: int_B={int_B} --> c_A={c_A} ({len(donor_path_BA)} rxns)",
              flush=True)

        # --- Step 4: Assemble final cross-feeders ---
        A_remaining_rxns = union_paths(path_dict_A, Core, exclude_core=c_A)
        cross_A = np.unique(np.concatenate([A_remaining_rxns, donor_path_BA])).astype(int)

        print(f"  CrossFeeder A: {len(cross_A)} rxns | CrossFeeder B: {len(cross_B)} rxns", flush=True)
        print(f"  A donates met {int_A} --> B uses for core {c_B}", flush=True)
        print(f"  B donates met {int_B} --> A uses for core {c_A}", flush=True)

        # --- Step 5: Verify viability of cross-feeders ---
        a_viable, a_missing = verify_autonomy(
            cross_A, rxnMat, prodMat, sumRxnVec,
            nutrientSet + [int_B], Currency, Core)
        b_viable, b_missing = verify_autonomy(
            cross_B, rxnMat, prodMat, sumRxnVec,
            nutrientSet + [int_A], Currency, Core)

        if not a_viable:
            print(f"  Viability check failed for cross_A. Missing cores: {a_missing}", flush=True)
            continue
        if not b_viable:
            print(f"  Viability check failed for cross_B. Missing cores: {b_missing}", flush=True)
            continue

        print(f"  Viability check passed.", flush=True)

        return {
            'auto_A':     full_rxns_A,
            'auto_B':     full_rxns_B,
            'cross_A':    cross_A,
            'cross_B':    cross_B,
            'A_donated':  int_A,
            'B_donated':  int_B,
            'A_ext_core': int(c_A),
            'B_ext_core': int(c_B),
            'pathway_AB': donor_pathway_AB,   # int_A --> c_B
            'pathway_BA': donor_path_BA,   # int_B --> c_A
        }

    print(f"\nERROR: Failed to construct cross-feeding pair after {max_attempts} attempts.",
          flush=True)
    return None


# ----------------------------------------------
# Test Run
# ----------------------------------------------

if __name__ == "__main__":
    import time
    from load_data import *
    from load_networks import load_minpaths

    # ── Config ──────────────────────────────────────────────────────────────
    PATHS_VERSION  = "1"
    USE_BYPRODUCTS = False
    MAX_ATTEMPTS   = 10
    # ────────────────────────────────────────────────────────────────────────

    all_paths = load_minpaths(f"paths_pv{PATHS_VERSION}")

    start = time.time()
    crossfeeders = build_crossfeeding_pair_from_paths(
        all_paths, rxnMat, prodMat, sumRxnVec,
        nutrientSet, Currency, Core,
        use_byproducts=USE_BYPRODUCTS,
        max_attempts=MAX_ATTEMPTS)
    end = time.time()
    print(f"\nTime taken: {end - start:.2f} seconds")

    if crossfeeders is not None:
        print("\nCross-feeding pair details:")
        print(f"    A needs metabolite {crossfeeders['B_donated']} from B to produce {crossfeeders['A_ext_core']}")
        print(f"    B needs metabolite {crossfeeders['A_donated']} from A to produce {crossfeeders['B_ext_core']}")
        print(f"    A --> B donor pathway has {len(crossfeeders['pathway_AB'])} reactions")
        print(f"    B --> A donor pathway has {len(crossfeeders['pathway_BA'])} reactions")
        print(f"    Net A: {len(crossfeeders['auto_A'])} --> {len(crossfeeders['cross_A'])} reactions")
        print(f"    Net B: {len(crossfeeders['auto_B'])} --> {len(crossfeeders['cross_B'])} reactions")
    else:
        print("\nFailed to construct cross-feeding pair.")
