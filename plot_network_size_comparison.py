import matplotlib.pyplot as plt
import numpy as np
from load_networks import *

'''
Plot the size distribution of autonomous and cross-feeding networks

Input format: dict with keys - type (auto/cross), subdir, file, label
'''

networks = [
    # {"type": "auto",  "subdir": "autonets_rs_av1",  "file": "P",       "label": "Reverse Scope Minimal Net"},
    # {"type": "auto",  "subdir": "autonets_mp_av1",  "file": "P_pv2",   "label": "Pathway Minimal Net"},
    # {"type": "auto",  "subdir": "autonets_mp_av1",  "file": "NP_pv2",  "label": "Pathway Non-Minimal Net"},
    {"type": "cross", "subdir": "crossnets_mp_cv1", "file": "int_P",   "label": "Minimal Cross-feeders"},
    {"type": "cross", "subdir": "crossnets_mp_cv1", "file": "byp_NP",  "label": "Non-Minimal Cross-feeders"},
]

all_sizes = []
labels = []
colors = ['#d0021b', '#4a90e2', '#7ed321', '#f5a623']

for d in networks:
    if d['type'] == 'auto':
        path = resolve_autonet_path(d['subdir'], d['file'])
    else:
        path = resolve_crossnet_path(d['subdir'], d['file'])
    label = d['label']
    with open(path, "rb") as f:
        results = pickle.load(f)
    if d['type'] == 'auto':
        sizes = np.array([len(net) for net in results])
    else:
        sizes = np.concatenate([
            np.array([len(p['cross_A']) for p in results]),
            np.array([len(p['cross_B']) for p in results]),
        ])
    all_sizes.append(sizes)
    labels.append(f"{label} $(\\mu={sizes.mean():.1f})$")
    print(f"{label}: N={len(results)}, min={sizes.min()}, max={sizes.max()}, mean={sizes.mean():.1f}")

fig, ax = plt.subplots(figsize=(6, 4))
for i, (sizes, label) in enumerate(zip(all_sizes, labels)):
    ax.hist(sizes, bins=30, alpha=0.7, weights=np.ones(len(sizes)) / len(sizes), 
            label=label, color=colors[i], histtype='stepfilled', edgecolor='k')
    ax.axvline(sizes.mean(), linestyle='--', linewidth=1, color=colors[i])
ax.set_xlabel("Network Size (# of reactions)")
ax.set_ylabel("Fraction of Networks")
ax.set_title("Network Size Distribution")
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()