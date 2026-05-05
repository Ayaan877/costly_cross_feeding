import pickle
import numpy as np
from directory_paths import *

'''
Load data for pre-generated pathways, networks and yields
'''

CORE_METS = [
    "C00009", "C00013", "C00022", "C00025",
    "C00041", "C00065", "C00097", "C00117"]


def load_minpaths(paths_subdir):
    """
    Load minimal pathways for all 8 core metabolites.
    paths_subdir : "paths_pv{pv}"  e.g. "paths_pv2"
    """
    all_targets = []
    for met in CORE_METS:
        path = resolve_paths_path(paths_subdir, met)
        with open(path, "rb") as f:
            results = pickle.load(f)
        all_targets.append(list(results['networks']))
    return all_targets


def load_autonets(autonet_subdir, autonet_file):
    """
    Load autonomous networks.

    autonet_subdir : "autonets_{source}_av{av}"
    autonet_file   : "{P|NP}_pv{pv}"  (mp)  or  "P"  (rs)
    """
    path = resolve_autonet_path(autonet_subdir, autonet_file)
    with open(path, "rb") as f:
        nets = pickle.load(f)
    print(f"Loaded {len(nets)} autonomous networks from {path}")
    return nets


def load_crossnets(crossnet_subdir, crossnet_file):
    """
    Load cross-feeding networks.

    crossnet_subdir : "crossnets_{source}_cv{cv}"
    crossnet_file   : "{byp|int}_{P|NP}"
    """
    path = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    with open(path, "rb") as f:
        nets = pickle.load(f)
    print(f"Loaded {len(nets)} cross-feeding networks from {path}")
    return nets


def load_yields(autonet_subdir, autonet_file, yield_mode,
                crossnet_subdir=None, crossnet_file=None):
    """
    Load networks and their yield data.

    autonet_subdir  : "autonets_{source}_av{av}"
    autonet_file    : "{P|NP}_pv{pv}"  (mp)  or  "P"  (rs)
    yield_mode      : "sbd" | "iter" | "alt" | "stoich"
    crossnet_subdir : "crossnets_{source}_cv{cv}"  (None → auto yield)
    crossnet_file   : "{byp|int}_{P|NP}"            (None → auto yield)

    Returns (net_sizes, E_yields, B_yields, via, nets).
    """
    typ = 'cross' if crossnet_subdir is not None else 'auto'
    net_path   = (resolve_crossnet_path(crossnet_subdir, crossnet_file)
                  if typ == 'cross'
                  else resolve_autonet_path(autonet_subdir, autonet_file))
    yield_path = resolve_yield_path(autonet_subdir, autonet_file, yield_mode,
                                    crossnet_subdir, crossnet_file)
    with open(net_path, "rb") as f:
        nets = pickle.load(f)
    with open(yield_path, "rb") as f:
        yields = pickle.load(f)

    if typ == "auto":
        if yield_mode == "stoich":
            via = np.ones(len(nets), dtype=bool)
            return np.array([len(n) for n in nets]), yields['fitness'], yields['biomass'], via, nets
        E_y, B_y, via = yields
        net_sizes = np.array([len(n) for n in nets])[via]
        return net_sizes, E_y[via], B_y[via], via, nets
    else:
        if yield_mode == "stoich":
            via  = np.ones(len(nets), dtype=bool)
            EA, EB = yields['fitness_A'], yields['fitness_B']
            BA, BB = yields['biomass_A'], yields['biomass_B']
        else:
            via  = yields['pair_viable']
            EA, EB = yields['E_A'][via], yields['E_B'][via]
            BA, BB = yields['B_A'][via], yields['B_B'][via]
        szA = np.array([len(p['cross_A']) for p in nets])[via]
        szB = np.array([len(p['cross_B']) for p in nets])[via]
        return (np.concatenate([szA, szB]),
                np.concatenate([EA, EB]),
                np.concatenate([BA, BB]),
                via, nets)


def load_merged_yields(crossnet_subdir, crossnet_file, yield_mode):
    """
    Load merged (union) network sizes and their yield data.

    crossnet_subdir : "crossnets_{source}_cv{cv}"
    crossnet_file   : "{byp|int}_{P|NP}"
    yield_mode      : "sbd" | "iter"

    Returns (net_sizes, E_yields, B_yields, via, nets).
    """
    net_path   = resolve_crossnet_path(crossnet_subdir, crossnet_file)
    yield_path = resolve_merged_yield_path(crossnet_subdir, crossnet_file, yield_mode)
    with open(net_path, "rb") as f:
        nets = pickle.load(f)
    with open(yield_path, "rb") as f:
        E_yields, B_yields, via = pickle.load(f)
    sizes = np.array([len(set(p['cross_A']) | set(p['cross_B'])) for p in nets])
    net_sizes = sizes[via]
    print(f"Loaded merged yields from {yield_path}")
    return net_sizes, E_yields[via], B_yields[via], via, nets
