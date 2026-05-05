'''
Sweep variants of toynet2 with R5 coupled to C (R5: P + a*C -> b*BA*),
optionally adjusting C output of R2.  Reports cross vs merged yield.
'''
import numpy as np
from calculate_crossNet_yield import splitByDemand_crossfeeding
from calculate_autoNet_yield import splitByDemand

# Metabolites:  X(0)  C(1)  P(2)  Q(3)  BA*(4)  AA*(5)
nutrientSet = [0]; Energy = []; Currency = []; Core = [4, 5]
Net1 = [0, 3]; Net2 = [1, 2, 4]; mergedNet = [0, 1, 2, 3, 4]
crossPair = {'cross_A': Net1, 'cross_B': Net2,
             'A_donated': 2, 'B_donated': 3}

def build(R2_C=2, R5_C=1, R5_BA=1, R5_P=1):
    S = np.array([
        [-1,  0,    3,   0,  3, 0],   # R1: X -> 3 BA* + 3 P
        [-1, R2_C,  0,   0,  0, 0],   # R2: X -> R2_C C
        [-1, -2,    0,   3,  0, 3],   # R3: X + 2C -> 3 AA* + 3 Q
        [ 0,  0,    0,  -1,  0, 1],   # R4: Q -> 1 AA*
        [ 0, -R5_C,-R5_P, 0, R5_BA,0],# R5: R5_P P + R5_C C -> R5_BA BA*
    ], dtype=float)
    return S

def run(S, label):
    rho = S.clip(max=0.0); pi = S.clip(min=0.0)
    rxnMat = (rho != 0).astype(int); prodMat = (pi != 0).astype(int)
    sumRxnVec = rxnMat.sum(axis=1)
    res = splitByDemand_crossfeeding(S, rxnMat, prodMat, sumRxnVec,
                                     rho, pi, nutrientSet, Energy,
                                     Currency, Core, crossPair)
    BA_cf, BB_cf = res['B_A'], res['B_B']
    via_cf = res['viable_A'] and res['viable_B']
    _, B_m, via_m = splitByDemand(S, rxnMat, prodMat, sumRxnVec,
                                  rho, pi, nutrientSet, Energy,
                                  Currency, Core, mergedNet)
    avg = (BA_cf + BB_cf) / 2
    ratio = avg / B_m if B_m else float('inf')
    flag = '  <-- OUTPERFORMS' if avg > B_m and via_cf and via_m else ''
    print(f"{label:35s} cf:({BA_cf:5.2f},{BB_cf:5.2f}) avg={avg:5.2f}  "
          f"merged={B_m:5.2f}  ratio={ratio:.2f}x  "
          f"viable cf={via_cf} merged={via_m}{flag}")

print("Baseline (no C coupling on R5):")
run(build(R2_C=2, R5_C=0, R5_BA=1, R5_P=1), "R5: P -> BA*")

print("\nCouple C to R5 (vary C demand / R2 output / R5 yield):")
for R2_C in [2, 3, 4, 5, 6]:
    for R5_C in [1, 2]:
        for R5_BA in [1, 2, 3]:
            lbl = f"R2:X->{R2_C}C  R5:P+{R5_C}C->{R5_BA}BA*"
            run(build(R2_C=R2_C, R5_C=R5_C, R5_BA=R5_BA), lbl)
