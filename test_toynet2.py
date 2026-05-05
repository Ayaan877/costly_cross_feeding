'''
TOY NETWORK 2 - OBLIGATE CROSS-FEEDING (intermediate-only exchange)

A "proper" obligate cross-feeder: neither organism exports core biomass
precursors. They exchange only non-core intermediate metabolites (P and Q).
Each organism still produces both BA* and AA* internally — but it can only
produce one of them efficiently from X alone, and needs the partner's
intermediate to make the other.

Network:
  Metabolites: X (nutrient), C, P, Q (intermediates), BA*, AA* (core)
  R1 (A):  X         →  3 BA* + 3 P          (efficient BA*; P is byproduct)
  R2 (B):  X         →  2 C                  (B's intermediate-only step)
  R3 (B):  X + 2 C   →  3 AA* + 3 Q          (efficient AA*; Q is byproduct)
  R4 (A):  Q         →  1 AA*                (A makes AA* from B's Q)
  R5 (B):  P         →  1 BA*                (B makes BA* from A's P)

Cross-feeding:
  Net A = {R1, R4}     donates P (intermediate) to Net B
  Net B = {R2, R3, R5} donates Q (intermediate) to Net A
  Exchanged metabolites P and Q are NON-core intermediates — no precursor
  exchange. Each organism still must internally synthesise both BA* and AA*.

Obligacy:
  Net A alone : R1 makes BA* but R4 has no Q source → no AA* → not viable.
  Net B alone : R3 makes AA* but R5 has no P source → no BA* → not viable.

Key nonlinearity (split-by-demand):
  In the merged network, X is split equally between R1, R2 and R3 (each
  demands 1 X), so R1 only gets 1/3 X → 1 BA* + 1 P, and R3 likewise yields
  1 AA* + 1 Q. The intermediate-conversion steps R4/R5 then give 1 more AA*
  and 1 more BA* → total biomass = 4.
  In cross-feeding, compartment A spends its full X on R1 → 3 BA* + 3 P,
  and compartment B's two X-using reactions split X 0.5/0.5, fully feeding
  R3 → 1.5 AA* + 1.5 Q. After exchange, each organism enjoys an efficient
  primary biomass route plus a "topped up" core met from the partner's
  intermediate, beating the merged yield.
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# ---------------------------------------------------------------------------
# Stoichiometric matrix  (rows = reactions, cols = metabolites)
# Metabolites:  X(0)  C(1)  P(2)  Q(3)  BA*(4)  AA*(5)
# Reactions:    R1(0) R2(1) R3(2) R4(3) R5(4)
# ---------------------------------------------------------------------------
stoich_matrix = np.array([
    [-1,  0,  3,  0,  3,  0],   # R1: X → 3 BA* + 3 P
    [-1,  2,  0,  0,  0,  0],   # R2: X → 2 C
    [-1, -2,  0,  3,  0,  3],   # R3: X + 2C → 3 AA* + 3 Q
    [ 0,  0,  0, -1,  0,  1],   # R4: Q → 1 AA*
    [ 0,  -1, -1,  0,  1,  0],   # R5: P + C → 1 BA*
], dtype=float)

rho         = stoich_matrix.clip(max=0.0)          # reactant-only stoich
pi          = stoich_matrix.clip(min=0.0)          # product-only stoich
rxnMat      = (rho != 0).astype(int)               # binary reactant matrix
prodMat     = (pi  != 0).astype(int)               # binary product matrix
sumRxnVec   = np.sum(rxnMat, axis=1)               # reactant count per rxn

nutrientSet = [0]        # X
Energy      = []
Currency    = []
Core        = [4, 5]     # BA* and AA*

met_names = ['X', 'C', 'P', 'Q', 'BA*', 'AA*']
rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5']

Net1       = [0, 3]       # R1, R4   (organism A)
Net2       = [1, 2, 4]    # R2, R3, R5 (organism B)
mergedNet  = [0, 1, 2, 3, 4]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 2,      # P (intermediate) flows A → B
    'B_donated': 3,      # Q (intermediate) flows B → A
}

print("=" * 50)
print("TOY NETWORK 2")
print("=" * 50)
print(f"Net A reactions : {[rxn_names[r] for r in Net1]}  →  donates P (intermediate)")
print(f"Net B reactions : {[rxn_names[r] for r in Net2]}  →  donates Q (intermediate)")
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

    ratio = (B1 + B2) / 2 / merged_B
    print(f"Cross-feeding / Merged  = {ratio:.2f}x")
    if (B1 + B2) / 2 > merged_B:
        print(">>> Cross-feeding OUTPERFORMS merged <<<")
    else:
        print("No outperformance.")
