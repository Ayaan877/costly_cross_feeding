'''
TOY NETWORK 7 — OBLIGATE CROSS-FEEDING WITH BRANCH-COMPETITION RELIEF
======================================================================
The merged network has TWO different reactions producing the same core
metabolite AA*, each using a different intermediate (P or Q). When
merged, both AA*-producers PLUS the BA*-producer PLUS both side-branch
reactions (R5,R6) all compete for the shared central intermediate I,
splitting flux 5 ways.

In each cross-feeding compartment, only HALF of those competing sinks
exist, so I is split fewer ways — yielding more biomass per unit X.
Each organism still needs the partner-derived metabolite to make AA*,
so each organism is obligately dependent on the other.

This avoids the "duplicate-reaction collapsed in merged" quirk because
the asymmetry is not in R4 (which is two distinct reactions R3 and R4
in merged), but in the broader sink-count for I.

Metabolites
-----------
  X    (0) — sole nutrient
  I    (1) — central intermediate (made from X)
  P    (2) — exchanged: A donates P to B
  Q    (3) — exchanged: B donates Q to A
  BA*  (4) — core biomass precursor 1
  AA*  (5) — core biomass precursor 2

Reactions
---------
  R1: X      → 2 I       (BOTH)   nutrient assimilation
  R2: I      → BA*       (BOTH)   makes core 1
  R3: I + Q  → AA*       (A only) A makes core 2 from Q (Q from B)
  R4: I + P  → AA*       (B only) B makes core 2 from P (P from A)
  R5: I      → P         (A only) A produces P (donates to B)
  R6: I      → Q         (B only) B produces Q (donates to A)

Standalone viability
--------------------
  Net A = {R1, R2, R3, R5}: makes I, BA*, P. R3 needs Q — NO Q producer.
        Cannot make AA*. NOT VIABLE.
  Net B = {R1, R2, R4, R6}: makes I, BA*, Q. R4 needs P — NO P producer.
        Cannot make AA*. NOT VIABLE.

When crossfeeding, each organism produces BOTH cores BA* and AA*.
Exchange: P (A→B) and Q (B→A) — both intermediates, never cores.

Why the cross-feeder wins (sink-count non-linearity)
----------------------------------------------------
MERGED (1 unit X total, all 6 reactions):
  Round 1: R1: X → 2 I.
  Round 2: I demanded by R2, R3, R4, R5, R6 — FIVE sinks each wanting 1.
           Each gets 2/5 = 0.4 I. R3, R4 wait (no Q, P yet).
           R2 → 0.4 BA*. R5 → 0.4 P. R6 → 0.4 Q.
  Round 3: R3 with 0.4 I + 0.4 Q → 0.4 AA*.
           R4 with 0.4 I + 0.4 P → 0.4 AA*.
  Merged total B = 0.4 + 0.4 + 0.4 = 1.2

CROSS-FEEDING (each org gets 1 unit X):
  Organism A {R1, R2, R3, R5} + exchange:
    Round 1: R1: X → 2 I.
    Round 2: I demanded by R2, R3, R5 — only THREE sinks.
             Each gets 2/3 ≈ 0.667 I. R3 waits for Q.
             R2 → 0.667 BA*. R5 → 0.667 P (exchanged to B).
    Round 3: Q arrives from B. R3 → 0.667 AA*.
  Organism A biomass ≈ 0.667 BA* + 0.667 AA* = 1.333.
  Organism B is symmetric.

  Cross-feeding average ≈ 1.333 vs merged 1.2  →  1.11× improvement.

This is robust to the "duplicate-reaction merging" issue because both
R3 and R4 are *distinct* reactions even in the merged net, and the
non-linearity comes purely from how many sinks share I.
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoich matrix (rows = reactions, cols = metabolites)
#  Metabolites:  X(0)  I(1)  P(2)  Q(3)  BA*(4)  AA*(5)
#  Reactions:    R1(0) R2(1) R3(2) R4(3) R5(4) R6(5)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  2,  0,  0,  0,  0],   # R1: X → 2I        (both)
    [ 0, -1,  0,  0,  1,  0],   # R2: I → BA*        (both)
    [ 0, -1,  0, -1,  0,  1],   # R3: I + Q → AA*    (A only)
    [ 0, -1, -1,  0,  0,  1],   # R4: I + P → AA*    (B only)
    [ 0, -1,  1,  0,  0,  0],   # R5: I → P          (A only)
    [ 0, -1,  0,  1,  0,  0],   # R6: I → Q          (B only)
], dtype=float)

rho       = stoich_matrix.clip(max=0.0)
pi        = stoich_matrix.clip(min=0.0)
rxnMat    = (rho != 0).astype(int)
prodMat   = (pi  != 0).astype(int)
sumRxnVec = np.sum(rxnMat, axis=1)

nutrientSet = [0]    # X
Energy      = []
Currency    = []
Core        = [4, 5] # BA* and AA*

rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
met_names = ['X', 'I', 'P', 'Q', 'BA*', 'AA*']

Net1      = [0, 1, 2, 4]   # Net A: {R1, R2, R3, R5} — needs Q from B
Net2      = [0, 1, 3, 5]   # Net B: {R1, R2, R4, R6} — needs P from A
mergedNet = [0, 1, 2, 3, 4, 5]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 2,    # P (intermediate) flows A → B
    'B_donated': 3,    # Q (intermediate) flows B → A
}

print("=" * 60)
print("TOY NETWORK 7 — branch-competition relief via cross-feeding")
print("=" * 60)
print("  R1: X → 2I        (both)")
print("  R2: I → BA*        (both)")
print("  R3: I + Q → AA*    (A only)  -- A's AA* route, needs Q from B")
print("  R4: I + P → AA*    (B only)  -- B's AA* route, needs P from A")
print("  R5: I → P          (A only)  -- A donates P")
print("  R6: I → Q          (B only)  -- B donates Q")
print()
print("Net A = {R1, R2, R3, R5}  receives Q from B, donates P")
print("Net B = {R1, R2, R4, R6}  receives P from A, donates Q")
print("Exchange: P and Q (intermediates — BA*/AA* never exchanged)")
print()

if __name__ == "__main__":
    # --- Standalone yields ---
    _, BA_solo, okA = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                    rho, pi, nutrientSet, Energy, Currency,
                                    Core, Net1)
    _, BB_solo, okB = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                    rho, pi, nutrientSet, Energy, Currency,
                                    Core, Net2)
    print(f"Net A alone              : B={BA_solo}  viable={okA}  (obligately dependent)")
    print(f"Net B alone              : B={BB_solo}  viable={okB}  (obligately dependent)")
    print()

    # --- Cross-feeding yields ---
    results = splitByDemand_crossfeeding(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, crossPair)

    B1    = results['B_A']
    B2    = results['B_B']
    stat1 = results['viable_A']
    stat2 = results['viable_B']

    print(f"Cross-feeding viability  : A={stat1}, B={stat2}")
    print(f"Cross-feeding yields     : Net A = {B1:.4f},  Net B = {B2:.4f}")
    print(f"  Average per organism   = {(B1 + B2) / 2:.4f}")
    print()

    # --- Merged network yield ---
    _, merged_B, merged_stat = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, mergedNet)

    print(f"Merged viability         : {merged_stat}")
    print(f"Merged yield             : {merged_B:.4f}")
    print()

    avg_cross = (B1 + B2) / 2
    if merged_B and merged_B > 0:
        ratio = avg_cross / merged_B
        print(f"Cross-feeding / Merged   = {ratio:.4f}x")
        if avg_cross > merged_B:
            print(">>> Cross-feeding OUTPERFORMS merged "
                  "(branch-competition relief) <<<")
        else:
            print("No outperformance.")
