'''
TOY NETWORK 3 - INTERMEDIATE EXCHANGE CROSS-FEEDING
Shows cross-feeding yield > merged yield, where only INTERMEDIATES are exchanged
(not core biomass precursors directly).

Network:
  Metabolites: X (nutrient), C (intermediate), P (intermediate), Q (intermediate),
               BA* and AA* (core biomass precursors — never exchanged)

  R1: X      → 3P        (Net A: efficient P producer, stoich 3)
  R2: X      → 2C        (Net B: C producer)
  R3: X + 2C → 3Q        (Net B: uses C + X to make Q)
  R4: P      → BA*       (shared: P → BA*)
  R5: Q      → AA*       (shared: Q → AA*)

Cross-feeding:
  Net A = {R1, R4, R5}  — specialises in making P; donates P to B; receives Q from B
  Net B = {R2, R3, R4, R5} — makes Q via the 2-step route; donates Q to A; receives P from A
  Exchange: A_donated = P (idx 2),  B_donated = Q (idx 3)
  Neither organism is viable without the partner (A has no Q-producing reaction;
  B has no P-producing reaction).

Key nonlinearity — same mechanism as toynet2, but one layer deeper:
  Merged: X split 3-ways among R1, R2, R3.  R1 only gets 1/3 X → 1P → 1 BA*.
  Cross:  A devotes full X to R1 → 3P; B devotes X only to R2+R3 (2-way split).
  Net A receives 0.75 Q from B → 0.75 AA*, plus 1.5 BA* of its own.
  Net B receives 1.5 P from A → 1.5 BA*, plus 0.75 AA* of its own.

Expected results:
  Merged yield       = 2.0  (1 BA* + 1 AA*)
  Cross-feeding avg  = 2.25 (1.5 BA* + 0.75 AA*  /  0.75 AA* + 1.5 BA*)
  Ratio              = 1.125x  — cross-feeding outperforms merged
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoichiometric matrix  (rows = reactions, cols = metabolites)
#
#  Metabolites:  X(0)  C(1)  P(2)  Q(3)  BA*(4)  AA*(5)
#  Reactions:    R1(0) R2(1) R3(2) R4(3) R5(4)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  0,  3,  0,  0,  0],   # R1: X → 3P
    [-1,  2,  0,  0,  0,  0],   # R2: X → 2C
    [-1, -2,  0,  3,  0,  0],   # R3: X + 2C → 3Q
    [ 0,  0, -1,  0,  1,  0],   # R4: P → BA*
    [ 0,  0,  0, -1,  0,  1],   # R5: Q → AA*
], dtype=float)

rho         = stoich_matrix.clip(max=0.0)
pi          = stoich_matrix.clip(min=0.0)
rxnMat      = (rho != 0).astype(int)
prodMat     = (pi  != 0).astype(int)
sumRxnVec   = np.sum(rxnMat, axis=1)

nutrientSet = [0]        # X
Energy      = []
Currency    = []
Core        = [4, 5]     # BA* and AA*

met_names = ['X', 'C', 'P', 'Q', 'BA*', 'AA*']
rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5']

# Net A specialises in P production; needs Q from B for AA*
Net1      = [0, 3, 4]        # R1, R4, R5
# Net B makes Q via the 2-step C route; needs P from A for BA*
Net2      = [1, 2, 3, 4]     # R2, R3, R4, R5
mergedNet = [0, 1, 2, 3, 4]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 2,      # P (intermediate) flows A → B
    'B_donated': 3,      # Q (intermediate) flows B → A
}

print("=" * 55)
print("TOY NETWORK 3 — intermediate exchange cross-feeding")
print("=" * 55)
print(f"Net A: {[rxn_names[r] for r in Net1]}   donates P (intermediate)")
print(f"Net B: {[rxn_names[r] for r in Net2]}  donates Q (intermediate)")
print(f"Exchange metabolites: P and Q  (NOT core mets BA*/AA*)")
print()

if __name__ == "__main__":
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
        print(">>> Cross-feeding OUTPERFORMS merged <<<")
    else:
        print("No outperformance.")
