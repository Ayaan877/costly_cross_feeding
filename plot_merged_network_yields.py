import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from directory_paths import *

'''
Compare the distribution of biomass yield differences between cross-feeding pairs
and their merged networks, for both byproduct and intermediate cross-feeding types. 
Highlights the fraction of pairs where the cross-feeding pair outperforms the merged 
network (positive delta) versus underperforms (negative delta).
'''

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
AUTONET_SUBDIR  = "autonets_mp_av1"
AUTONET_FILE    = "P_pv2"
CROSSNET_SUBDIR = "crossnets_mp_cv1"
CROSSNET_FILES  = ("byp_P", "int_P")
YIELD_MODE      = "sbd"

COLOR_POS  = "C0"   # highlight: crossfeeder > merged
COLOR_NEG  = "C1"   # crossfeeder <= merged
ALPHA      = 0.85
BINS       = 30

CTYPE_LABELS = {"byp": "byproduct", "int": "intermediate"}

# ---------------------------------------------------------------------------
# Resolve comparison file path (written by outperformer_vs_merged.py)
# ---------------------------------------------------------------------------
source, av, pruning, _ = parse_autonet_spec(AUTONET_SUBDIR, AUTONET_FILE)
_, cv, _, _            = parse_crossnet_spec(CROSSNET_SUBDIR, CROSSNET_FILES[0])

comp_path = YIELDS_DIR / f"yields_{source}_{YIELD_MODE}" / \
            f"comparison_vs_merged_{pruning}_av{av}_cv{cv}.pkl"

print(f"Loading: {comp_path}")
with open(comp_path, "rb") as f:
    results = pickle.load(f)

# ---------------------------------------------------------------------------
# Build shared bin edges across both ctypes so axes are comparable
# delta_pair = (B_A + B_B) - 2 * merged_B  (one value per valid pair)
# ---------------------------------------------------------------------------
all_deltas = []
for crossnet_file in CROSSNET_FILES:
    _, _, ctype, _ = parse_crossnet_spec(CROSSNET_SUBDIR, crossnet_file)
    res = results[ctype]
    vm  = res['valid_mask']
    pair_delta = res['delta_B'][vm]
    all_deltas.append(pair_delta)

all_deltas_cat = np.concatenate(all_deltas)
lo, hi  = np.nanmin(all_deltas_cat), np.nanmax(all_deltas_cat)
# Force 0 to be a bin edge by splitting BINS proportionally between sides
n_neg = max(1, round(BINS * (-lo) / (hi - lo)))
n_pos = max(1, BINS - n_neg)
bin_edges = np.concatenate([
    np.linspace(lo, 0, n_neg + 1),
    np.linspace(0, hi, n_pos + 1)[1:],
])

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, len(CROSSNET_FILES), figsize=(12, 4.5), sharey=False)
fig.suptitle(
    r"$\Delta$ Biomass yield: $(B_A + B_B)_{\mathrm{cross}}$ $-$ $B_{\mathrm{merged}}$",
    fontsize=13)

for ax, crossnet_file in zip(axes, CROSSNET_FILES):
    _, _, ctype, _ = parse_crossnet_spec(CROSSNET_SUBDIR, crossnet_file)
    res  = results[ctype]
    vm   = res['valid_mask']
    n_valid = int(np.sum(vm))

    delta = res['delta_B'][vm]
    n_pos = int(np.sum(delta > 0))

    counts, edges, patches = ax.hist(
        delta, bins=bin_edges,
        weights=np.ones(len(delta)) / len(delta),
        color=COLOR_NEG, alpha=ALPHA, edgecolor='none')

    # Highlight bars whose left edge is >= 0 (outperformer region)
    for patch, left in zip(patches, edges[:-1]):
        if left >= 0:
            patch.set_facecolor(COLOR_POS)

    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')

    ax.legend(handles=[
        Patch(color=COLOR_POS, alpha=ALPHA, label=f"Outperformers ({n_pos/n_valid:.1%})"),
        Patch(color=COLOR_NEG, alpha=ALPHA, label=f"Non-outperformers ({(n_valid - n_pos)/n_valid:.1%})"),
    ], fontsize=10)

    ax.set_title(
        f"{CTYPE_LABELS[ctype].capitalize()} Cross-Feeders\n"
        f"$N_{{pairs}}={n_valid}$",
        fontsize=10)
    ax.set_xlabel(r"$\Delta$ Biomass Yield", fontsize=11)
    ax.set_ylabel("Fraction of Pairs", fontsize=11)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
