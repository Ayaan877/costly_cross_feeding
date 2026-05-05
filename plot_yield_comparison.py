import matplotlib.pyplot as plt
import numpy as np
from load_networks import *

'''
Plot comparison of yields for autonomous, cross-feeding and merged networks

Input format: dict with keys - type (auto/cross/cross_merged), subdir, file, yield_mode, label
  type "auto"         : autonet subdir + file
  type "cross"        : crossnet subdir + file (requires autonet subdir + file for paired autonets)
  type "cross_merged" : crossnet subdir + file, merged yields
'''
networks = [
    {"type":"auto", "subdir":"autonets_mp_av1", "file":"P_pv2", "yield_mode":"sbd", "label":"Autonomous Networks"},
    {"type":"cross", "subdir":"crossnets_mp_cv1", "file":"byp_P", "yield_mode":"sbd", "label":"Byproduct Crossfeeders"},
    {"type":"cross", "subdir":"crossnets_mp_cv1", "file":"int_P", "yield_mode":"sbd", "label":"Intermediate Crossfeeders"},
    # {"type":"cross_merged", "subdir":"crossnets_mp_cv1", "file":"int_P", "yield_mode":"sbd", "label":"Merged Intermediate Cross-feeders"},
    # {"type":"cross_merged", "subdir":"crossnets_mp_cv1", "file":"byp_P", "yield_mode":"sbd", "label":"Merged Byproduct Cross-feeders"},
]

all_data = []
for d in networks:
    label      = d['label']
    yield_mode = d['yield_mode']
    if d['type'] == 'cross_merged':
        net_sizes, viable_E, viable_B, via, nets = load_merged_yields(
            d['subdir'], d['file'], yield_mode)
    elif d['type'] == 'cross':
        net_sizes, viable_E, viable_B, via, nets = load_yields(
            None, None, yield_mode, d['subdir'], d['file'])
    else:  # auto
        net_sizes, viable_E, viable_B, via, nets = load_yields(
            d['subdir'], d['file'], yield_mode, None, None)
    all_data.append((f"{label}", net_sizes, viable_E, viable_B))
    print(f"{label}: N={len(nets)}, viable={len(viable_E)}, "
          f"size mean={net_sizes.mean():.1f}, max E={viable_E.max():.3f}, max B={viable_B.max():.3f}")

colors = ['#d0021b', '#4a90e2', '#7ed321', '#f5a623']
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for i, (label, net_sizes, viable_E, viable_B) in enumerate(all_data):
    axes[0].hist(net_sizes, bins=20, alpha=0.5, weights=np.ones(len(net_sizes))/len(net_sizes), 
                 color=colors[i], label=label, histtype='stepfilled', edgecolor='k')
    axes[0].axvline(net_sizes.min(), color=colors[i], linestyle='--', linewidth=1)

    axes[1].hist(viable_E,  bins=20, alpha=0.5, weights=np.ones(len(viable_E))/len(viable_E), 
                 color=colors[i], label=label, histtype='stepfilled', edgecolor='k')
    axes[1].axvline(viable_E.max(),  color=colors[i], linestyle='--', linewidth=1)

    axes[2].hist(viable_B,  bins=20, alpha=0.5, weights=np.ones(len(viable_B))/len(viable_B), 
                 color=colors[i], label=label, histtype='stepfilled', edgecolor='k')
    axes[2].axvline(viable_B.max(),  color=colors[i], linestyle='--', linewidth=1)

axes[0].set_xlabel("Network Size (# reactions)"); axes[0].set_ylabel("Fraction of Networks"); axes[0].legend(fontsize=10, loc='upper right')
axes[1].set_xlabel("Energy Yield");               axes[1].set_ylabel("Fraction of Networks"); axes[1].legend(fontsize=10, loc='upper right')
axes[2].set_xlabel("Biomass Yield");              axes[2].set_ylabel("Fraction of Networks"); axes[2].legend(fontsize=10, loc='upper right')
fig.suptitle("Network Yield - Byproduct vs Intermediate Cross-feeders")
plt.tight_layout()
plt.show()