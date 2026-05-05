'''
TOY NETWORK 1 - CROSS-FEEDING WITH SHARED INTERMEDIATES

Two organisms that each cover one efficient route to a core metabolite but
depend on intermediates from the partner to complete the other.

Network:
  Metabolites: X (nutrient), A, B, BAA, AA (intermediates),
               A_star, BA_star, AA_star (core)
  R1 (A):  X            →  3A + B
  R2 (A):  A + B         →  BA_star
  R3 (B):  A             →  A_star
  R4 (B):  BAA           →  AA_star
  R5 (B):  X             →  A + BAA
  R6 (B):  BAA           →  B + AA
  R7 (A):  2A            →  AA_star
  R8 (B):  B + A_star    →  BA_star

Cross-feeding:
  Net A = {R1, R2, R7}        donates A_star  to Net B
  Net B = {R3, R4, R5, R6, R8} donates AA     to Net A
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoichiometric matrix  (rows = reactions, cols = metabolites)
# Metabolites:  X(0)  A(1)  B(2)  BAA(3)  AA(4)  A_star(5)  BA_star(6)  AA_star(7)
# Reactions:    R1(0) R2(1) R3(2) R4(3)   R5(4)  R6(5)      R7(6)       R8(7)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  3,  1,  0,  0,  0,  0,  0],   # R1: X → 3A + B
    [ 0, -1, -1,  0,  0,  0,  1,  0],   # R2: A + B → BA_star
    [ 0, -1,  0,  0,  0,  1,  0,  0],   # R3: A → A_star
    [ 0,  0,  0, -1,  0,  0,  0,  1],   # R4: BAA → AA_star
    [-1,  1,  0,  1,  0,  0,  0,  0],   # R5: X → A + BAA
    [ 0,  0,  1, -1,  1,  0,  0,  0],   # R6: BAA → B + AA
    [ 0, -2,  0,  0,  0,  0,  0,  1],   # R7: 2A → AA_star
    [ 0,  0, -1,  0,  0, -1,  1,  0],   # R8: B + A_star → BA_star
], dtype=float)

rho         = stoich_matrix.clip(max=0.0)
pi          = stoich_matrix.clip(min=0.0)
rxnMat      = (rho != 0).astype(int)
prodMat     = (pi  != 0).astype(int)
sumRxnVec   = np.sum(rxnMat, axis=1)

nutrientSet = [0]        # X
Energy      = []
Currency    = []
Core        = [6, 7]     # BA_star, AA_star

met_names = ['X', 'A', 'B', 'BAA', 'AA', 'A_star', 'BA_star', 'AA_star']
rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']

Net1      = [0, 1, 6]          # R1, R2, R7  (organism A)
Net2      = [2, 3, 4, 5, 7]    # R3, R4, R5, R6, R8 (organism B)
mergedNet = list(range(8))

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 5,    # A_star flows A → B
    'B_donated': 4,    # AA flows B → A
}

print("=" * 50)
print("TOY NETWORK 1")
print("=" * 50)
print(f"Net A reactions : {[rxn_names[r] for r in Net1]}  →  donates A_star")
print(f"Net B reactions : {[rxn_names[r] for r in Net2]}  →  donates AA")
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

    print(f"Cross-feeding viability : {stat1 and stat2}")
    print(f"Cross-feeding yields    : Net A = {B1:.2f},  Net B = {B2:.2f}")
    print(f"  Average per nutrient  = {(B1 + B2) / 2:.2f}")
    print()

    # --- Merged network yield ---
    _, merged_B, merged_stat = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, mergedNet)

    print(f"Merged viability        : {merged_stat}")
    print(f"Merged yield            : {merged_B:.2f}")
    print()

    ratio = (B1 + B2) / 2 / merged_B if merged_B else float('inf')
    print(f"Cross-feeding / Merged  = {ratio:.2f}x")
    if (B1 + B2) / 2 > merged_B:
        print(">>> Cross-feeding OUTPERFORMS merged <<<")
    else:
        print("No outperformance.")
