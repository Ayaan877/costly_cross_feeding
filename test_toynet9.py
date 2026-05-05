'''
TOY NETWORK 9 - INTERMEDIATE-LIMITED + OBLIGATE + OUTPERFORMER

Combines three properties:
  (a) The merged network's biomass-producing reactions are limited by an
      intermediate metabolite (I), not the nutrient.
  (b) The cross-feeding pair is OBLIGATE — neither organism is viable alone
      (each makes only one core precursor; each needs the partner's core).
  (c) Cross-feeding average yield strictly exceeds merged yield (nonlinear
      "splitting metabolic work" advantage).

Network:
  Metabolites:  X(0) nutrient,  I(1) intermediate,  BA*(2), AA*(3) core
  R1:  X       →  4 I        (shared intermediate producer)
  R2:  I       →  3 BA*      (cheap, I-light BA* branch)
  R3:  4 I     →  2 AA*      (I-heavy AA* branch)

Per-compartment specialisation:
  Net A = {R1, R2}   – produces BA* only.  Donates BA* → B.
  Net B = {R1, R3}   – produces AA* only.  Donates AA* → A.
  Each compartment lacks one core precursor → both obligately depend on
  partner exchange to reach viability.

Merged-network split-by-demand trace:
  Round 1: R1 sole X-user  →  X-share = 1, limRct = X.  Produces 4 I.
  Round 2: R2 demands 1 I, R3 demands 4 I  →  total I demand = 5.
           R2 I-share = (1·4)/5 = 0.8;   limRct(R2) = I (intermediate)
                BA* = 0.8 · 3 = 2.4
           R3 I-share = (4·4)/5 = 3.2;   limRct(R3) = I (intermediate)
                AA* = (3.2/4) · 2 = 1.6
  Merged biomass total = 2.4 + 1.6 = 4.0
  Both biomass-producing reactions are limited by I, not X.

Cross-feed split-by-demand trace:
  Compartment A: R1 X-share = 1 → 4 IA.  R2 sole I-demander → uses all 4 IA.
                 limRct(R2) = I (still I-limited, but no longer split).
                 BA*A = 4 · 3 = 12.
  Compartment B: R1 X-share = 1 → 4 IB.  R3 sole I-demander → uses all 4 IB.
                 limRct(R3) = I, AA*B = (4/4)·2 = 2.
  Exchange: BA* split 50/50 (A→B), AA* split 50/50 (B→A).
            Each organism sees:  6 BA* + 1 AA* = 7 biomass.
  Cross-feed average = 7.0  vs  merged 4.0   →  1.75× advantage.

Bottleneck identity:
  Merged    : R2 and R3 are both I-limited; I is split 1:4 between them.
  Cross-fed : R2 (in A) and R3 (in B) are still I-limited, but each receives
              the FULL I supply of its compartment — the intermediate-level
              competition is dissolved by specialisation.
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# Metabolites:  X(0)  I(1)  BA*(2)  AA*(3)
# Reactions:    R1(0) R2(1) R3(2)
stoich_matrix = np.array([
    [-1,  4,  0,  0],   # R1: X → 4 I
    [ 0, -1,  3,  0],   # R2: I → 3 BA*
    [ 0, -4,  0,  2],   # R3: 4 I → 2 AA*
], dtype=float)

rho       = stoich_matrix.clip(max=0.0)
pi        = stoich_matrix.clip(min=0.0)
rxnMat    = (rho != 0).astype(int)
prodMat   = (pi  != 0).astype(int)
sumRxnVec = np.sum(rxnMat, axis=1)

nutrientSet = [0]
Energy      = []
Currency    = []
Core        = [2, 3]

met_names = ['X', 'I', 'BA*', 'AA*']
rxn_names = ['R1', 'R2', 'R3']

Net1      = [0, 1]      # A: R1, R2  (BA* maker; needs AA* from B)
Net2      = [0, 2]      # B: R1, R3  (AA* maker; needs BA* from A)
mergedNet = [0, 1, 2]

crossPair = {
    'cross_A':   Net1,
    'cross_B':   Net2,
    'A_donated': 2,     # BA* A → B
    'B_donated': 3,     # AA* B → A
}

print("=" * 64)
print("TOY NETWORK 9 — intermediate-limited, obligate, outperformer")
print("=" * 64)
print(f"Net A : {[rxn_names[r] for r in Net1]}  (makes BA*, needs AA*)")
print(f"Net B : {[rxn_names[r] for r in Net2]}  (makes AA*, needs BA*)")
print()

if __name__ == "__main__":
    # --- Obligacy check: each net alone ---
    _, B_A_alone, stat_A_alone = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, Net1)
    _, B_B_alone, stat_B_alone = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, Net2)
    print(f"Net A alone viable : {stat_A_alone}   (must be False for obligacy)")
    print(f"Net B alone viable : {stat_B_alone}   (must be False for obligacy)")
    print()

    # --- Merged network ---
    _, merged_B, merged_stat = splitByDemand(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, mergedNet)
    print(f"Merged viability  : {merged_stat}")
    print(f"Merged yield      : {merged_B:.4f}")
    print()

    # --- Cross-feeding ---
    results = splitByDemand_crossfeeding(
        stoich_matrix, rxnMat, prodMat, sumRxnVec,
        rho, pi, nutrientSet, Energy, Currency,
        Core, crossPair)
    B1, B2 = results['B_A'], results['B_B']
    print(f"Cross-feed viable : {results['viable_A'] and results['viable_B']}")
    print(f"Cross-feed yields : Net A = {B1:.4f}, Net B = {B2:.4f}")
    print(f"  Average         = {(B1 + B2) / 2:.4f}")
    print()

    ratio = (B1 + B2) / 2 / merged_B if merged_B else float('nan')
    print(f"Cross-feed / Merged = {ratio:.3f}x")
    if (B1 + B2) / 2 > merged_B:
        print(">>> Cross-feeding OUTPERFORMS merged <<<")
    print()
    print("Limiting-reactant identity:")
    print("  Merged    : R2 and R3 both limited by intermediate I "
          "(I split 1:4 between them).")
    print("  Cross-fed : R2 in A and R3 in B each see the full I supply"
          " of their compartment — intermediate-level competition removed.")
