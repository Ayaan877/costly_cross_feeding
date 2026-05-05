'''
TOY NETWORK 6 — OBLIGATE CROSS-FEEDING THAT PERFORMS WORSE THAN MERGED
=======================================================================
Both organisms are obligately dependent (neither is viable alone).
The merged network OUTPERFORMS cross-feeding because compartmentalisation
forces the exchanged intermediate to be split between self-use and donation,
whereas the merged network uses it all directly — relieving a bottleneck
that cross-feeding itself creates.

Metabolites
-----------
  X   (0) — sole nutrient
  C   (1) — intermediate  [exchanged: B donates C to A]
  P   (2) — intermediate  [exchanged: A donates P to B]
  BA* (3) — core biomass precursor 1
  AA* (4) — core biomass precursor 2

Reactions
---------
  R1: X  →  C        (Net B only: nutrient into intermediate C)
  R2: C  →  P        (Net A only: converts C into P)
  R3: P  →  BA*      (both nets: P → biomass precursor 1)
  R4: C  →  AA*      (both nets: C → biomass precursor 2)

Viability
---------
  Net A = {R2, R3, R4}: needs C — no C producer. NOT VIABLE ALONE.
  Net B = {R1, R3, R4}: needs P for R3 — no P producer. NOT VIABLE ALONE.

Exchange: C (B → A),  P (A → B)  — both intermediates, never core mets.

Why merged OUTPERFORMS cross-feeding (flux bottleneck in the exchange)
----------------------------------------------------------------------
MERGED — single pool, all reactions share metabolites directly:

  Round 1: X available. Only R1 fires (R2 needs C; R3 needs P; R4 needs C).
    X=1 → 1C

  Round 2: C available. R2 and R4 both fire, splitting C equally (demand 1:1):
    R2 gets 0.5C → 0.5P
    R4 gets 0.5C → 0.5 AA*

  Round 3: P available. R3 fires:
    0.5P → 0.5 BA*

  Merged B = 0.5 BA* + 0.5 AA* = 1.0

CROSS-FEEDING — C_B is produced in compartment B:

  Round 1: B's R1 fires: X_B → C_B = 1

  Round 2: C_B is split THREE ways — B's R4, B's exchange (→C_A), and nothing else.
    Wait: C_B is demanded by B's R4 AND the exchange reaction (B→A).
    Both have demand 1 → split 50/50:
      B's R4 gets 0.5C_B → 0.5 AA*_B
      Exchange fires: 0.5C_B → 0.5C_A

  Round 3: C_A = 0.5 available in A. A's R2 and A's R4 split C_A equally:
    A's R2 gets 0.25C_A → 0.25P_A
    A's R4 gets 0.25C_A → 0.25 AA*_A

  Round 4: P_A = 0.25 split between A's R3 and exchange (P: A→B):
    A's R3 gets 0.125P → 0.125 BA*_A
    Exchange: 0.125P_A → 0.125P_B

  Round 5: B's R3 fires: 0.125P_B → 0.125 BA*_B

  Cross-feeding yields:
    B_A = 0.125 BA* + 0.25 AA* = 0.375
    B_B = 0.125 BA* + 0.5  AA* = 0.625
    Average = 0.5 < 1.0 merged

THE BOTTLENECK: In cross-feeding, the exchange reaction competes with B's R4
for C_B (Round 2: 3-way split becomes 2-way split but with exchange taking half).
This halves the C available to R4 in B, and also halves the C sent to A.
The intermediate then loses another factor of 2 at each downstream split.
Each exchange step is a "tax" on the intermediate — the more steps in the
chain, the more flux is lost to these splits.

In merged, there is no exchange step: C goes directly to both R2 and R4 in
the same round, without paying the exchange tax. The merged pool allows the
intermediate to be shared without the overhead of an explicit transfer reaction.

Expected outputs
----------------
  Net A alone   : not viable
  Net B alone   : not viable
  Merged yield  : 1.0
  Cross-feeding : avg = 0.5  (A=0.375, B=0.625)
  Merged wins by 2×
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoich matrix  (rows = reactions, cols = metabolites)
#  Metabolites:  X(0)  C(1)  P(2)  BA*(3)  AA*(4)
#  Reactions:    R1(0) R2(1) R3(2) R4(3)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  1,  0,  0,  0],   # R1: X → C        (Net B only)
    [ 0, -1,  1,  0,  0],   # R2: C → P        (Net A only)
    [ 0,  0, -1,  1,  0],   # R3: P → BA*       (both nets)
    [ 0, -1,  0,  0,  1],   # R4: C → AA*       (both nets)
], dtype=float)

rho       = stoich_matrix.clip(max=0.0)
pi        = stoich_matrix.clip(min=0.0)
rxnMat    = (rho != 0).astype(int)
prodMat   = (pi  != 0).astype(int)
sumRxnVec = np.sum(rxnMat, axis=1)

nutrientSet = [0]    # X
Energy      = []
Currency    = []
Core        = [3, 4] # BA* and AA*

Net1      = [1, 2, 3]      # Net A: {R2, R3, R4}  — needs C from B
Net2      = [0, 2, 3]      # Net B: {R1, R3, R4}  — needs P from A
mergedNet = [0, 1, 2, 3]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 2,    # P (intermediate) flows A → B
    'B_donated': 1,    # C (intermediate) flows B → A
}

print("=" * 60)
print("TOY NETWORK 6 — obligate cross-feeding WORSE than merged")
print("=" * 60)
print("  R1: X → C        (Net B only)")
print("  R2: C → P        (Net A only)")
print("  R3: P → BA*       (both nets)")
print("  R4: C → AA*       (both nets)")
print()
print("Net A = {R2, R3, R4}  receives C from B, donates P")
print("Net B = {R1, R3, R4}  receives P from A, donates C")
print("Exchange: C and P (intermediates — BA*/AA* never exchanged)")
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
    print(f"  Average per nutrient   = {(B1 + B2) / 2:.4f}")
    print()

    # --- Merged network yield ---
    _, merged_B, merged_stat = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, mergedNet)

    print(f"Merged viability         : {merged_stat}")
    print(f"Merged yield             : {merged_B:.4f}")
    print()

    ratio = merged_B / ((B1 + B2) / 2)
    print(f"Merged / Cross-feeding   = {ratio:.4f}x")
    if merged_B > (B1 + B2) / 2:
        print(">>> Merged OUTPERFORMS cross-feeding (exchange-splitting bottleneck) <<<")
    else:
        print("No merged advantage.")
