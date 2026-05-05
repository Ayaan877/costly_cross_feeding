"""
Plots the cumulative frequency of pairwise reaction overlap (Jaccard index).

One panel per crossnet type (byp / int), shown side by side.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from directory_paths import *

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CROSSNET_SUBDIR = "crossnets_mp_cv1"
CROSSNET_FILES  = ("byp_P", "int_P")   # both crossnet variants

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def jaccard(a, b):
    sa, sb = set(a), set(b)
    union  = sa | sb
    return len(sa & sb) / len(union) if union else 0.0


def ecdf(values):
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------
CTYPE_LABELS = {"byp": "Byproduct", "int": "Intermediate"}

COLORS = {
    "autonomous" : "C0",
    "crossfeed"  : "C1",
}
LABELS = {
    "autonomous" : "Autonomous network pairs",
    "crossfeed"  : "Cross-feeding network pairs",
}
LW = 1.8

# ---------------------------------------------------------------------------
# Build figure  (one panel per crossnet type, side by side)
# ---------------------------------------------------------------------------
n_types = len(CROSSNET_FILES)
fig, axes = plt.subplots(1, n_types, figsize=(6 * n_types, 5), sharey=True)
if n_types == 1:
    axes = [axes]

for ax, crossnet_file in zip(axes, CROSSNET_FILES):
    _, _, ctype, _ = parse_crossnet_spec(CROSSNET_SUBDIR, crossnet_file)

    crossnet_path = resolve_crossnet_path(CROSSNET_SUBDIR, crossnet_file)
    print(f"\n[{ctype}] CrossNets : {crossnet_path}")
    with open(crossnet_path, "rb") as f:
        CrossNets = pickle.load(f)

    indices = range(len(CrossNets))

    jac_auto  = np.array([jaccard(CrossNets[i]['auto_A'],  CrossNets[i]['auto_B'])  for i in indices])
    jac_cross = np.array([jaccard(CrossNets[i]['cross_A'], CrossNets[i]['cross_B']) for i in indices])

    print(f"  Total pairs          : {len(CrossNets)}")
    print(f"  Auto  median Jaccard : {np.median(jac_auto):.3f}")
    print(f"  Cross median Jaccard : {np.median(jac_cross):.3f}")

    series = [
        (jac_auto,  "autonomous"),
        (jac_cross, "crossfeed"),
    ]
    panel_title = f"{CTYPE_LABELS.get(ctype, ctype)} cross-feeders"

    # ------------------------------------------------------------------
    # Plot CDFs
    # ------------------------------------------------------------------
    for jac, key in series:
        if len(jac) == 0:
            continue
        x, y = ecdf(jac)
        ax.step(x, y, where="post", color=COLORS[key], label=LABELS[key], linewidth=LW)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Reaction overlap (Jaccard index)", fontsize=11)
    ax.set_title(panel_title, fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=9, loc="lower right")

axes[0].set_ylabel("Cumulative frequency", fontsize=11)

fig.suptitle("Pairwise reaction overlap", fontsize=12)
plt.tight_layout()
plt.show()
