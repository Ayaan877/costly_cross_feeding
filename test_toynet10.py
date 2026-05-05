'''
TOY NETWORK 10 — INTERMEDIATE-EXCHANGE OBLIGATE CROSS-FEEDER
                  (intermediate-limited merged, cross-feed outperforms)

Properties:
  (1) Exchange is INTERMEDIATE-only.  Donated mets are P, Q (non-core).
      No core biomass precursor is ever transferred between organisms.
  (2) Pair is OBLIGATE — neither organism is viable in isolation:
        Net A alone : R4 needs P (B-only) → no AA* → not viable.
        Net B alone : R5 needs Q (A-only) → no BA* → not viable.
  (3) The merged net's biomass-producing reactions R2, R3 are limited by
      the intermediate I (not by X), and cross-feeding strictly OUTPERFORMS
      merged because the intermediate-level demand competition is reduced.

Network:
  Mets:  X(0)  I(1)  P(2)  Q(3)  BA*(4)  AA*(5)
  R1:  X     → 4 I              (shared X→I producer)
  R2:  I     → 3 BA*            (A only — A's main biomass branch)
  R3:  4 I   → 2 AA*            (B only — B's main biomass branch, I-heavy)
  R4:  P     → 1 AA*            (A only — A's AA* via partner's P)
  R5:  Q     → 1 BA*            (B only — B's BA* via partner's Q)
  R6:  I     → 1 Q              (A only — A produces Q intermediate)
  R7:  I     → 1 P              (B only — B produces P intermediate)

Networks:
  Net A = {R1, R2, R4, R6}   — donates Q,  receives P
  Net B = {R1, R3, R5, R7}   — donates P,  receives Q
  Merged = union {R1..R7}

Merged trace (split-by-demand):
  Round 1: only R1 has its reactant (X). R1 fires, X-share=1, limRct=X,
           produces 4 I.
  Round 2: I demand from R2(1) + R3(4) + R6(1) + R7(1) = 7.
           Shares (of 4 I produced):
             R2 = 4·1/7 = 4/7  → BA* = (4/7)·3       = 12/7
             R3 = 4·4/7 = 16/7 → AA* = (16/7)/4·2    =  8/7   (limRct = I)
             R6 = 4/7          → Q produced          =  4/7
             R7 = 4/7          → P produced          =  4/7
  Round 3: R4 sees P=4/7 → AA* = 4/7;   R5 sees Q=4/7 → BA* = 4/7
  Merged biomass = 12/7 + 8/7 + 4/7 + 4/7 = 28/7 = 4.000
  >>> R2 and R3 are both limited by intermediate I (I split 4-ways). <<<

Cross-feed trace:
  Compartment A {R1, R2, R4, R6}:
    R1 → 4 IA. I demand in A = R2(1) + R6(1) = 2.
    R2 share = 4·1/2 = 2  → BA*A = 6
    R6 share = 2          → 2 QA produced
    R4 awaits PA from exchange.
  Compartment B {R1, R3, R5, R7}:
    R1 → 4 IB. I demand in B = R3(4) + R7(1) = 5.
    R3 share = 4·4/5 = 3.2  → AA*B = (3.2/4)·2 = 1.6   (limRct = I)
    R7 share = 4·1/5 = 0.8  → 0.8 PB produced
    R5 awaits QB from exchange.
  Exchange (full transfer to partner because no in-compartment demand):
    QA = 2   →  R5 in B uses 1 Q step → BA*B = 2
    PB = 0.8 →  R4 in A uses 0.8 P → AA*A = 0.8
  Net A biomass = 6 + 0.8 = 6.8
  Net B biomass = 1.6 + 2  = 3.6
  Average = 5.2   →  Cross-feed / Merged = 5.2 / 4.0 = 1.30×
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# Mets:  X(0)  I(1)  P(2)  Q(3)  BA*(4)  AA*(5)
# Rxns:  R1(0) R2(1) R3(2) R4(3) R5(4) R6(5) R7(6)
stoich_matrix = np.array([
    [-1,  4,  0,  0,  0,  0],   # R1: X → 4 I
    [ 0, -1,  0,  0,  3,  0],   # R2: I → 3 BA*
    [ 0, -4,  0,  0,  0,  2],   # R3: 4 I → 2 AA*
    [ 0,  0, -1,  0,  0,  1],   # R4: P → 1 AA*
    [ 0,  0,  0, -1,  1,  0],   # R5: Q → 1 BA*
    [ 0, -1,  0,  1,  0,  0],   # R6: I → 1 Q
    [ 0, -1,  1,  0,  0,  0],   # R7: I → 1 P
], dtype=float)

rho       = stoich_matrix.clip(max=0.0)
pi        = stoich_matrix.clip(min=0.0)
rxnMat    = (rho != 0).astype(int)
prodMat   = (pi  != 0).astype(int)
sumRxnVec = np.sum(rxnMat, axis=1)

nutrientSet = [0]
Energy      = []
Currency    = []
Core        = [4, 5]

met_names = ['X', 'I', 'P', 'Q', 'BA*', 'AA*']
rxn_names = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']

Net1      = [0, 1, 3, 5]   # A: R1, R2, R4, R6  → donates Q, receives P
Net2      = [0, 2, 4, 6]   # B: R1, R3, R5, R7  → donates P, receives Q
mergedNet = [0, 1, 2, 3, 4, 5, 6]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 3,        # Q (intermediate) flows A → B
    'B_donated': 2,        # P (intermediate) flows B → A
}

print("=" * 68)
print("TOY NETWORK 10 — intermediate-exchange obligate cross-feeder")
print("=" * 68)
print(f"Net A : {[rxn_names[r] for r in Net1]}  donates Q (intermediate)")
print(f"Net B : {[rxn_names[r] for r in Net2]}  donates P (intermediate)")
print()

if __name__ == "__main__":
    # --- Obligacy check ---
    _, _, A_alone = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                  rho, pi, nutrientSet, Energy, Currency,
                                  Core, Net1)
    _, _, B_alone = splitByDemand(stoich_matrix, rxnMat, prodMat, sumRxnVec,
                                  rho, pi, nutrientSet, Energy, Currency,
                                  Core, Net2)
    print(f"Net A alone viable : {A_alone}   (must be False)")
    print(f"Net B alone viable : {B_alone}   (must be False)")
    print()

    # --- Merged ---
    _, merged_B, merged_stat = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, mergedNet)
    print(f"Merged viability  : {merged_stat}")
    print(f"Merged yield      : {merged_B:.4f}")
    print()

    # --- Cross-feed ---
    results = splitByDemand_crossfeeding(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, crossPair)
    B1, B2 = results['B_A'], results['B_B']
    print(f"Cross-feed viable : {results['viable_A'] and results['viable_B']}")
    print(f"Cross-feed yields : Net A = {B1:.4f}, Net B = {B2:.4f}")
    print(f"  Average         = {(B1 + B2) / 2:.4f}")
    print(f"  flux Q (A→B)    = {results['flux_A_to_B']:.4f}")
    print(f"  flux P (B→A)    = {results['flux_B_to_A']:.4f}")
    print()

    ratio = (B1 + B2) / 2 / merged_B if merged_B else float('nan')
    print(f"Cross-feed / Merged = {ratio:.3f}x")
    if (B1 + B2) / 2 > merged_B:
        print(">>> Cross-feeding OUTPERFORMS merged <<<")
    print()
    print("Limiting-reactant identity:")
    print("  Merged    : R2, R3 both limRct = I  (I split 4-ways across "
          "R2/R3/R6/R7).")
    print("  Cross-fed : in A, I split only between R2 and R6 (2-way);"
          " in B, I split between R3 and R7. Intermediate-level demand"
          " competition reduced → outperformance.")
