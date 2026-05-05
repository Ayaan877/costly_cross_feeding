import matplotlib.pyplot as plt
import numpy as np
from load_networks import *

'''
Plot pathway discovery rate, size distribution, and similarity across targets
'''

PATHS_VERSION = "2"
files = ["C00009", "C00013", "C00022", "C00025", "C00041", "C00065", "C00097", "C00117"]

fig, axes = plt.subplots(2, 4, figsize=(12, 6))
for ax, file in zip(axes.flat, files):
    with open(resolve_paths_path(f"paths_pv{PATHS_VERSION}", file), "rb") as f:
        results = pickle.load(f)

    paths = results['networks']
    unique_counts = np.array(results['unique_counts'])
    attempts = np.array(results['attempts'])*8

    ax.plot(attempts, unique_counts,'o-')
    ax.set_title(f"{file} ({len(paths)} paths)")
    ax.set_xlabel("# of Attempts")
    ax.set_ylabel("Unique Pathways")
    ax.grid()
fig.suptitle(f"Unique Pathway Discovery Rate (pv{PATHS_VERSION})")
plt.tight_layout()
plt.show()

fig2, axes2 = plt.subplots(2, 4, figsize=(12, 6))
for ax, file in zip(axes2.flat, files):
    with open(resolve_paths_path(f"paths_pv{PATHS_VERSION}", file), "rb") as f:
        results = pickle.load(f)

    paths = results['networks']
    path_sizes = np.array([len(path) for path in paths])

    ax.hist(path_sizes, bins=20, alpha=0.7, weights=np.ones(len(path_sizes)) / len(path_sizes))
    ax.axvline(path_sizes.mean(), color='red', linestyle='--', label=f"Mean: {path_sizes.mean():.1f}")
    ax.set_title(f"{file} ({len(path_sizes)} paths)")
    ax.set_xlabel("# of Reactions")
    ax.set_ylabel("Fraction of Pathways")
    ax.legend()
fig2.suptitle(f"Pathway Size Distribution (pv{PATHS_VERSION})")
plt.tight_layout()
plt.show()

fig3, axes3 = plt.subplots(2, 4, figsize=(12, 6))
for ax, file in zip(axes3.flat, files):
    with open(resolve_paths_path(f"paths_pv{PATHS_VERSION}", file), "rb") as f:
        results = pickle.load(f)

    paths = results['networks']
    n_paths = len(paths)

    # Count how often each reaction appears across all pathways for this target
    rxn_counts = {}
    for path in paths:
        for rxn in path:
            rxn_counts[rxn] = rxn_counts.get(rxn, 0) + 1

    # Fraction of pathways containing each reaction
    rxn_fracs = np.array(list(rxn_counts.values())) / n_paths

    ax.hist(rxn_fracs, bins=20, alpha=0.7)
    ax.set_title(f"{file} ({len(rxn_counts)} rxns)")
    ax.set_xlabel("Proportion of Pathways")
    ax.set_ylabel("# of Reactions Shared")
    ax.set_yscale('log')
    ax.grid(which='both', alpha=0.5)
fig3.suptitle(f"Reaction Sharing Across Pathways (pv{PATHS_VERSION})")
plt.tight_layout()
plt.show()