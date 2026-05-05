import matplotlib.pyplot as plt
import numpy as np
from load_networks import *

'''
Check Network Yield

Input format: dict with keys - type (auto/cross/cross_merged), subdir, file, yield_mode, label
    type "auto"         : autonet subdir + file
    type "cross"        : crossnet subdir + file (requires autonet subdir + file for paired autonets)
    type "cross_merged" : crossnet subdir + file, merged yields 
'''

dataset = {"type": "auto", "subdir": "autonets_mp_av1", "file": "NP_pv2", "yield_mode": "sbd", "label": "Autonomous Non-Minimal Network Yields"}

if dataset['type'] == 'cross_merged':
    net_sizes, viable_E, viable_B, viability, nets = load_merged_yields(
        dataset['subdir'], dataset['file'], dataset['yield_mode'])
elif dataset['type'] == 'cross':
    net_sizes, viable_E, viable_B, viability, nets = load_yields(
        None, None, dataset['yield_mode'], dataset['subdir'], dataset['file'])
else:  # auto
    net_sizes, viable_E, viable_B, viability, nets = load_yields(
        dataset['subdir'], dataset['file'], dataset['yield_mode'], None, None)
label = dataset['label']
print(f"{label}: N={len(nets)}, viable={len(viable_E)}, "
      f"size mean={net_sizes.mean():.1f}, max E={viable_E.max():.3f}, max B={viable_B.max():.3f}")

viable_indices = np.arange(len(net_sizes))
min_size_idx = np.argmin(net_sizes)
max_E_idx    = np.argmax(viable_E)
max_B_idx    = np.argmax(viable_B)
colors = ['#d0021b', '#4a90e2', '#7ed321', '#f5a623']

fig, ax = plt.subplots(1, 3, figsize=(10, 4))
ax[0].hist(net_sizes, bins=20, histtype='stepfilled', alpha=0.8, facecolor=colors[0], weights=np.ones(len(net_sizes)) / len(net_sizes), edgecolor='black')
ax[0].axvline(net_sizes[min_size_idx], color=colors[0], lw=2, label="Smallest Network")
ax[0].axvline(net_sizes[max_E_idx],    color=colors[1], lw=2, label="Highest Energy Yield Network")
ax[0].axvline(net_sizes[max_B_idx],    color=colors[2], lw=2, label="Highest Biomass Yield Network")
ax[0].set_xlabel("Network Size (# reactions)"); ax[0].set_ylabel("Fraction of Networks")

ax[1].hist(viable_E, bins=20, histtype='stepfilled', alpha=0.8, facecolor=colors[1], weights=np.ones(len(viable_E)) / len(viable_E), edgecolor='black')
ax[1].axvline(viable_E[min_size_idx], color=colors[0], lw=2, label="Smallest Network")
ax[1].axvline(viable_E[max_E_idx],    color=colors[1], lw=2, label="Highest Energy Yield Network")
ax[1].axvline(viable_E[max_B_idx],    color=colors[2], lw=2, label="Highest Biomass Yield Network")
ax[1].set_xlabel("Energy Yield"); ax[1].set_ylabel("Fraction of Networks")

ax[2].hist(viable_B, bins=20, histtype='stepfilled', alpha=0.8, facecolor=colors[2], weights=np.ones(len(viable_B)) / len(viable_B), edgecolor='black')
ax[2].axvline(viable_B[min_size_idx], color=colors[0], lw=2, label="Smallest Network")
ax[2].axvline(viable_B[max_E_idx],    color=colors[1], lw=2, label="Highest Energy Yield Network")
ax[2].axvline(viable_B[max_B_idx],    color=colors[2], lw=2, label="Highest Biomass Yield Network")
ax[2].set_xlabel("Biomass Yield"); ax[2].set_ylabel("Fraction of Networks")
handles, labels_leg = ax[0].get_legend_handles_labels()
fig.legend(handles, labels_leg, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.93), fontsize=10)
fig.suptitle(f"{label} ($N={len(viable_E)}$)", y=0.98, fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.93))
plt.show()
