'''
TOY NETWORK 4 — LIMITING-REACTANT SWITCH MECHANISM
=====================================================
This example shows that cross-feeding outperforms the merged network specifically
because MERGING CHANGES WHICH REACTANT IS LIMITING in a two-substrate reaction,
degrading the biomass yield.

Metabolites
-----------
  X   (0) — sole nutrient
  P   (1) — intermediate (exchanged: A donates P to B)
  Q   (2) — intermediate (exchanged: B donates Q to A)
  BA* (3) — core biomass precursor 1
  AA* (4) — core biomass precursor 2

Reactions
---------
  R1: X      →  4P          (Net A only:  efficient P-producing route, stoich 4)
  R2: X      →  P           (Net B only:  inefficient P-producing route, stoich 1)
  R3: X      →  2Q          (both nets:   Q producer)
  R4: 2P + Q →  3 BA*       (both nets:   P-heavy two-substrate synthesis)
  R5: P + 2Q →  3 AA*       (both nets:   Q-heavy two-substrate synthesis)

Network assignments
-------------------
  Net A = {R1, R3, R4, R5}  — uses efficient R1; P-rich; donates surplus P
  Net B = {R2, R3, R4, R5}  — uses inefficient R2; P-scarce; donates Q

Both organisms are VIABLE ALONE (each has R3 for Q, R4+R5 for core mets).
Cross-feeding lets B receive cheap P from A, removing P as the bottleneck in B.

Exchange: P (A→B),  Q (B→A)  — BOTH are intermediates, never core mets.

Limiting-reactant switch — worked example
------------------------------------------
Net B ALONE:  X split 2-ways (R2 and R3, equal demand):
  R2 gets 1/2 X → 1/2 P
  R3 gets 1/2 X → 1 Q

  P (1/2) and Q (1) are split between R4 (demands 2P:1Q) and R5 (demands 1P:2Q):
    total P-demand = 3,  total Q-demand = 3
    R4 P-share = 2/3 * 1/2 = 1/3,  R4 Q-share = 1/3 * 1 = 1/3
    R5 P-share = 1/3 * 1/2 = 1/6,  R5 Q-share = 2/3 * 1 = 2/3

    R4 (2P:1Q):  need 1/3 * (1/2) = 1/6 Q if P is limit, have 1/3 Q  → P is limiting
      fires at P-scale 1/3  → BA* = (1/3) * (3/2) = 0.5
    R5 (1P:2Q):  need 1/6 * 2 = 1/3 Q if P is limit, have 2/3 Q  → P is limiting
      fires at P-scale 1/6  → AA* = (1/6) * (3/1) = 0.5
    Net B alone B = 0.5 + 0.5 = 1.0

Net A ALONE:  X split 2-ways (R1 and R3):
  R1 gets 1/2 X → 2P    (stoich 4, so 4 * 1/2 = 2)
  R3 gets 1/2 X → 1Q

  P (2) and Q (1): R4 P-share = 4/3, Q-share = 1/3; R5 P-share = 2/3, Q-share = 2/3

    R4: need 4/3 * (1/2) = 2/3 Q if P is limit, have 1/3 Q  → Q is limiting
      fires at Q-scale 1/3  → BA* = (1/3) * (3/1) = 1.0
    R5: need 2/3 * 2 = 4/3 Q if P is limit, have 2/3 Q  → Q is limiting
      fires at Q-scale 2/3  → AA* = (2/3) * (3/2) = 1.0
    Net A alone B = 1.0 + 1.0 = 2.0

MERGED net = {R1, R2, R3, R4, R5}:  X split 3-ways (R1, R2, R3):
  R1 gets 1/3 X → 4/3 P,  R2 gets 1/3 X → 1/3 P,  R3 gets 1/3 X → 2/3 Q
  Total P = 4/3 + 1/3 = 5/3,  Q = 2/3

  R4 P-share = 2/3 * 5/3 = 10/9,  R4 Q-share = 1/3 * 2/3 = 2/9
  R5 P-share = 1/3 * 5/3 = 5/9,   R5 Q-share = 2/3 * 2/3 = 4/9

    R4: need 10/9 * 1/2 = 5/9 Q if P-limit, have 2/9  → Q is limiting!
      fires at Q-scale 2/9  → BA* = (2/9) * 3 = 2/3
    R5: need 5/9 * 2 = 10/9 Q if P-limit, have 4/9  → Q is limiting!
      fires at Q-scale 4/9  → AA* = (4/9) * (3/2) = 2/3
    Merged B = 2/3 + 2/3 = 4/3 ≈ 1.333

  THE SWITCH: Net B alone has P as limiting reactant in R4 and R5.
  Merged has Q as limiting reactant in R4 and R5.
  Adding R1 into the pool floods P but does NOT proportionally increase Q,
  so the limiting reactant SWITCHES from P to Q, and the richer Q units
  are now what gets wasted — giving LOWER yield than Net A alone.

Cross-feeding average = (B_A + B_B) / 2 = (2.0 + 1.0) / 2 = 1.5 > 1.333 merged.
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoich matrix  (rows = reactions, cols = metabolites)
#  Metabolites:  X(0)  P(1)  Q(2)  BA*(3)  AA*(4)
#  Reactions:    R1(0) R2(1) R3(2) R4(3)  R5(4)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  4,  0,  0,  0],   # R1: X → 4P    (Net A: efficient P)
    [-1,  1,  0,  0,  0],   # R2: X → P     (Net B: inefficient P)
    [-1,  0,  2,  0,  0],   # R3: X → 2Q   (both nets)
    [ 0, -2, -1,  3,  0],   # R4: 2P+Q → 3BA*
    [ 0, -1, -2,  0,  3],   # R5: P+2Q → 3AA*
], dtype=float)

rho       = stoich_matrix.clip(max=0.0)
pi        = stoich_matrix.clip(min=0.0)
rxnMat    = (rho != 0).astype(int)
prodMat   = (pi  != 0).astype(int)
sumRxnVec = np.sum(rxnMat, axis=1)

nutrientSet = [0]
Energy      = []
Currency    = []
Core        = [3, 4]   # BA* and AA*

met_names = ['X', 'P', 'Q', 'BA*', 'AA*']
rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5']

Net1      = [0, 2, 3, 4]   # R1, R3, R4, R5  (Net A — efficient P route)
Net2      = [1, 2, 3, 4]   # R2, R3, R4, R5  (Net B — inefficient P route)
mergedNet = [0, 1, 2, 3, 4]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 1,    # P flows A → B
    'B_donated': 2,    # Q flows B → A
}

print("=" * 60)
print("TOY NETWORK 4 — limiting-reactant switch mechanism")
print("=" * 60)
print("  R1: X → 4P          (Net A only — efficient P producer)")
print("  R2: X → P           (Net B only — inefficient P producer)")
print("  R3: X → 2Q          (both nets)")
print("  R4: 2P + Q → 3 BA*  (both nets — P-heavy)")
print("  R5: P + 2Q → 3 AA*  (both nets — Q-heavy)")
print("Net A = {R1,R3,R4,R5}  donates surplus P")
print("Net B = {R2,R3,R4,R5}  receives P from A, donates Q")
print("Exchange: P and Q  (intermediates — BA*/AA* never exchanged)")
print()

if __name__ == "__main__":
    # --- Standalone yields for reference ---
    _, BA_solo, okA = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                    rho, pi, nutrientSet, Energy, Currency,
                                    Core, Net1)
    _, BB_solo, okB = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                    rho, pi, nutrientSet, Energy, Currency,
                                    Core, Net2)
    print(f"Net A alone              : B={BA_solo:.4f}  viable={okA}")
    print(f"Net B alone              : B={BB_solo:.4f}  viable={okB}")
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

    ratio = (B1 + B2) / 2 / merged_B
    print(f"Cross-feeding / Merged   = {ratio:.4f}x")
    if (B1 + B2) / 2 > merged_B:
        print(">>> Cross-feeding OUTPERFORMS merged (limiting-reactant switch) <<<")
    else:
        print("No outperformance.")
