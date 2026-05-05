"""
Central path resolver for all data files.

Naming convention
-----------------
Paths      data/paths/paths_pv{pv}/paths_{core_id}.pkl
AutoNets   data/networks/autonets_{source}_av{av}/autonets_{P|NP}[_pv{pv}].pkl
CrossNets  data/networks/crossnets_{source}_cv{cv}/crossnets_{byp|int}_{P|NP}.pkl
Yields     data/yields/yields_{source}_{mode}/yields_{auto|cross}[_{byp|int}]_{P|NP}_{av|cv}{v}.pkl

Subdirectory + file strings
------------------------------------
autonet_subdir  : "autonets_{source}_av{av}"     e.g. "autonets_mp_av2"
autonet_file    : "{P|NP}_pv{pv}"  (mp)  or "P" (rs)
crossnet_subdir : "crossnets_{source}_cv{cv}"    e.g. "crossnets_mp_cv1"
crossnet_file   : "{byp|int}_{P|NP}"             e.g. "byp_NP"
paths_subdir    : "paths_pv{pv}"                 e.g. "paths_pv1"
"""

from pathlib import Path

NETWORKS_DIR = Path("data/networks")
YIELDS_DIR   = Path("data/yields")
PATHS_DIR    = Path("data/paths")


# ── Spec parsers ─────────────────────────────────────────────────────────────

def parse_autonet_spec(subdir, file_detail):
    """
    Parse autonet subdirectory and file-detail strings.

    subdir      : "autonets_{source}_av{av}"
    file_detail : "{P|NP}_pv{pv}"  (mp)  or  "P"  (rs)

    Returns (source, av, pruning, pv); pv is None for rs.
    """
    parts   = subdir.split("_")         # ["autonets", source, "av{av}"]
    source  = parts[1]
    av      = parts[2][2:]              # strip "av"
    fparts  = file_detail.split("_", 1)
    pruning = fparts[0]                 # P or NP
    pv      = fparts[1][2:] if len(fparts) > 1 else None   # strip "pv"
    return source, av, pruning, pv


def parse_crossnet_spec(subdir, file_detail):
    """
    Parse crossnet subdirectory and file-detail strings.

    subdir      : "crossnets_{source}_cv{cv}"
    file_detail : "{byp|int}_{P|NP}"

    Returns (source, cv, cross_type, pruning).
    """
    parts      = subdir.split("_")     # ["crossnets", source, "cv{cv}"]
    source     = parts[1]
    cv         = parts[2][2:]          # strip "cv"
    cross_type, pruning = file_detail.split("_", 1)
    return source, cv, cross_type, pruning


# ── Path resolvers ────────────────────────────────────────────────────────────

def resolve_paths_path(paths_subdir, core_id):
    """
    paths_subdir : "paths_pv{pv}"   e.g. "paths_pv1"
    core_id      : KEGG metabolite ID  e.g. "C00022"
    """
    return PATHS_DIR / paths_subdir / f"paths_{core_id}.pkl"


def resolve_autonet_path(autonet_subdir, autonet_file):
    """
    autonet_subdir : "autonets_{source}_av{av}"
    autonet_file   : "{P|NP}_pv{pv}"  (mp)  or  "P"  (rs)
    """
    _, _, pruning, pv = parse_autonet_spec(autonet_subdir, autonet_file)
    fname = f"autonets_{pruning}_pv{pv}.pkl" if pv else f"autonets_{pruning}.pkl"
    return NETWORKS_DIR / autonet_subdir / fname


def resolve_crossnet_path(crossnet_subdir, crossnet_file):
    """
    crossnet_subdir : "crossnets_{source}_cv{cv}"
    crossnet_file   : "{byp|int}_{P|NP}"
    """
    _, _, cross_type, pruning = parse_crossnet_spec(crossnet_subdir, crossnet_file)
    return NETWORKS_DIR / crossnet_subdir / f"crossnets_{cross_type}_{pruning}.pkl"


def resolve_yield_path(autonet_subdir, autonet_file, yield_mode,
                       crossnet_subdir=None, crossnet_file=None):
    """
    Resolve the path for a yield pickle file.

    autonet_subdir  : "autonets_{source}_av{av}"
    autonet_file    : "{P|NP}_pv{pv}"  or  "P"
    yield_mode      : "sbd" | "iter" | "alt" | "stoich"
    crossnet_subdir : "crossnets_{source}_cv{cv}"  (None → auto yield)
    crossnet_file   : "{byp|int}_{P|NP}"            (None → auto yield)
    """
    if crossnet_subdir is None:
        source, av, pruning, _ = parse_autonet_spec(autonet_subdir, autonet_file)
        subdir_name = f"yields_{source}_{yield_mode}"
        fname = f"yields_auto_{pruning}_av{av}.pkl"
    else:
        source_cross, cv, cross_type, cross_pruning = parse_crossnet_spec(
            crossnet_subdir, crossnet_file)
        if autonet_subdir is not None:
            source, _, _, _ = parse_autonet_spec(autonet_subdir, autonet_file)
        else:
            source = source_cross
        subdir_name = f"yields_{source}_{yield_mode}"
        fname = f"yields_cross_{cross_type}_{cross_pruning}_cv{cv}.pkl"
    return YIELDS_DIR / subdir_name / fname


def resolve_merged_yield_path(crossnet_subdir, crossnet_file, yield_mode):
    """
    Resolve the path for a merged-network yield pickle file.

    The merged network is the union of both organisms' reactions in a cross-
    feeding pair, treated as a single autonomous network for yield calculation.

    crossnet_subdir : "crossnets_{source}_cv{cv}"
    crossnet_file   : "{byp|int}_{P|NP}"
    yield_mode      : "sbd" | "iter"
    """
    source, cv, cross_type, cross_pruning = parse_crossnet_spec(
        crossnet_subdir, crossnet_file)
    subdir_name = f"yields_{source}_{yield_mode}"
    fname = f"yields_merged_{cross_type}_{cross_pruning}_cv{cv}.pkl"
    return YIELDS_DIR / subdir_name / fname