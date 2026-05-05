import pickle
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx

"""
Visualise a single autonomous metabolic network as an undirected metabolite graph.

Nodes are metabolites, edges connect any reactant-product pair sharing a
reaction in the network. Currency metabolites are omitted from the graph
and listed instead in a "Growth medium" legend box, mirroring the figure
style of the thesis.
"""

# Make every text element bold by default
plt.rcParams.update({
    "font.weight":       "bold",
    "axes.titleweight":  "bold",
    "axes.labelweight":  "bold",
    "figure.titleweight":"bold",
})

from load_data import *

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
AUTONET_PATH = Path("data/networks/autonets_mp_av1/autonets_P_pv2.pkl")
SELECT       = "smallest"     # "smallest" | int index
OUT_PATH     = Path("plots/autonet_smallest_mp_pv2.png")

TITLE        = "A representative autonomous metabolic network"

# Colours
COLOR_SEED   = "#7ed321"   # green
COLOR_INTER  = "#f5a623"   # orange / amber
COLOR_CORE   = "#d0021b"   # red
EDGE_COLOR   = "#000000"   # black

LABEL_NUTRIENTS = True
LABEL_CORES     = True
LABEL_INTERMEDIATES = False   # set True to label everything

SIZE_BASE   = 60     # node area for degree 1
SIZE_STEP   = 90     # additional area per extra degree

# ---------------------------------------------------------------------------
# Load network
# ---------------------------------------------------------------------------
with open(AUTONET_PATH, "rb") as f:
    autonets = pickle.load(f)

if SELECT == "smallest":
    sizes = np.array([len(n) for n in autonets])
    idx   = int(np.argmin(sizes))
else:
    idx   = int(SELECT)

net = np.asarray(autonets[idx])
print(f"Loaded {len(autonets)} networks from {AUTONET_PATH}")
print(f"Plotting network #{idx} with {len(net)} reactions")

# ---------------------------------------------------------------------------
# Build metabolite graph (currency metabolites excluded)
# ---------------------------------------------------------------------------
currency_set = set(Currency)
core_set     = set(Core)
nutrient_set = set(nutrientSet)

G = nx.Graph()

for r in net:
    reactants = [m for m in np.nonzero(rxnMat[r])[0] if m not in currency_set]
    products  = [m for m in np.nonzero(prodMat[r])[0] if m not in currency_set]
    for m in reactants + products:
        G.add_node(int(m))
    for a in reactants:
        for b in products:
            if a == b:
                continue
            G.add_edge(int(a), int(b))

print(f"Graph: {G.number_of_nodes()} metabolite nodes, "
      f"{G.number_of_edges()} edges")

# ---------------------------------------------------------------------------
# Node attributes
# ---------------------------------------------------------------------------
def node_color(m):
    if m in nutrient_set: return COLOR_SEED
    if m in core_set:     return COLOR_CORE
    return COLOR_INTER

def met_name(m):
    kegg = inv_met_map.get(m, f"idx{m}")
    name = cpd_string_dict.get(kegg, kegg)
    # take only the first synonym, trim
    name = name.split(";")[0].strip()
    if len(name) > 22:
        name = name[:20] + "..."
    return name

nodes      = list(G.nodes())
colors     = [node_color(m) for m in nodes]
degrees    = dict(G.degree())
node_sizes = [SIZE_BASE + SIZE_STEP * (degrees[m] - 1) for m in nodes]

# Labels
labels = {}
for m in nodes:
    if m in nutrient_set and LABEL_NUTRIENTS:
        labels[m] = met_name(m)
    elif m in core_set and LABEL_CORES:
        labels[m] = met_name(m)
    elif LABEL_INTERMEDIATES:
        labels[m] = met_name(m)

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
pos = nx.spring_layout(G, k=1.2 / np.sqrt(max(G.number_of_nodes(), 1)),
                       iterations=400, seed=7)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 8))

nx.draw_networkx_edges(G, pos, edge_color=EDGE_COLOR, width=0.7,
                       alpha=0.7, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors,
                       node_size=node_sizes, linewidths=0.4,
                       edgecolors="black", ax=ax)

# Offset labels above their node, scaled by node radius so they clear the circle.
# matplotlib node_size is the marker AREA in points^2, so radius (in points) = sqrt(area)/2.
# Convert points -> data units using the current axis transform.
fig.canvas.draw()
trans = ax.transData.inverted()

def _radius_data_units(area_pts2):
    # marker radius in points, then convert points -> pixels -> data units (y-axis)
    r_pts = np.sqrt(area_pts2) / 2.0
    r_px  = r_pts * fig.dpi / 72.0
    # Use the y-axis scale: convert a vertical pixel offset to data delta
    (_, y0) = trans.transform((0, 0))
    (_, y1) = trans.transform((0, r_px))
    return abs(y1 - y0)

size_lookup = dict(zip(nodes, node_sizes))
LABEL_PAD_PTS = 4   # extra gap between circle edge and label baseline (points)
pad_data = _radius_data_units(LABEL_PAD_PTS ** 2 * 4)  # ~ LABEL_PAD_PTS in data units

label_pos = {}
for m, (x, y) in pos.items():
    r = _radius_data_units(size_lookup.get(m, SIZE_BASE))
    label_pos[m] = (x, y + r + pad_data)

text_artists = nx.draw_networkx_labels(G, label_pos, labels=labels, font_size=10,
                                       font_color="black", font_weight="bold",
                                       ax=ax,
                                       verticalalignment='bottom')

# ---------------------------------------------------------------------------
# Push labels away from any node circle they overlap.
#
# For each label we compute its bounding box (display coords), check overlap
# with every node disc, and shift the label radially outward (relative to its
# anchor node) until it is clear of all nodes (or we hit MAX_ITERS).
# ---------------------------------------------------------------------------
fig.canvas.draw()  # ensure text bboxes are valid
renderer = fig.canvas.get_renderer()  # type: ignore[attr-defined]

# Pre-compute node centres + radii in display (pixel) coords
node_disp = {m: ax.transData.transform(pos[m]) for m in nodes}
node_rad_px = {}
for m in nodes:
    r_pts = np.sqrt(size_lookup[m]) / 2.0
    node_rad_px[m] = r_pts * fig.dpi / 72.0

EXTRA_PAD_PX = 3   # gap between circle edge and label bbox
LABEL_PAD_PX = 6   # gap between two labels
MAX_ITERS    = 200
STEP_PX      = 4

def _bbox_disp(text):
    return text.get_window_extent(renderer=renderer)

def _bbox_overlaps_any_node(bbox, anchor_m):
    # Closest point on bbox to circle centre; if distance < radius -> overlap
    for m in nodes:
        cx, cy = node_disp[m]
        r = node_rad_px[m] + EXTRA_PAD_PX
        nx_ = max(bbox.xmin, min(cx, bbox.xmax))
        ny_ = max(bbox.ymin, min(cy, bbox.ymax))
        dx = cx - nx_
        dy = cy - ny_
        if dx * dx + dy * dy < r * r:
            return m
    return None

def _bbox_overlaps_other_label(bbox, self_key, others):
    # `others` is dict {key: bbox}
    for k, b in others.items():
        if k == self_key:
            continue
        # Inflate by LABEL_PAD_PX
        if (bbox.xmin < b.xmax + LABEL_PAD_PX and
            bbox.xmax + LABEL_PAD_PX > b.xmin and
            bbox.ymin < b.ymax + LABEL_PAD_PX and
            bbox.ymax + LABEL_PAD_PX > b.ymin):
            # return centre of conflicting bbox so we can push away from it
            return ((b.xmin + b.xmax) / 2.0, (b.ymin + b.ymax) / 2.0)
    return None

inv_trans = ax.transData.inverted()

# Iterate a few global passes so that pushing one label can be re-checked
# against everything else (label-vs-label is order-dependent).
GLOBAL_PASSES = 4
for _pass in range(GLOBAL_PASSES):
    # snapshot current bboxes for label-label checks
    fig.canvas.draw()
    label_bboxes = {k: _bbox_disp(t) for k, t in text_artists.items()}
    moved_any = False

    for m, txt in text_artists.items():
        ax_x, ax_y = node_disp[m]
        lx_disp, ly_disp = ax.transData.transform(txt.get_position())
        # Default push direction: away from anchor node
        dxv = lx_disp - ax_x
        dyv = ly_disp - ax_y
        norm = np.hypot(dxv, dyv)
        if norm < 1e-6:
            dxv, dyv, norm = 0.0, 1.0, 1.0
        ux, uy = dxv / norm, dyv / norm

        for _ in range(MAX_ITERS):
            bbox = _bbox_disp(txt)
            hit_node  = _bbox_overlaps_any_node(bbox, m)
            hit_label = _bbox_overlaps_other_label(bbox, m, label_bboxes)
            if hit_node is None and hit_label is None:
                break

            # Choose push direction: away from the offending node or label
            if hit_node is not None and hit_node != m:
                hx, hy = node_disp[hit_node]
                cx = (bbox.xmin + bbox.xmax) / 2.0
                cy = (bbox.ymin + bbox.ymax) / 2.0
                dxv2, dyv2 = cx - hx, cy - hy
                n2 = np.hypot(dxv2, dyv2)
                if n2 > 1e-6:
                    ux, uy = dxv2 / n2, dyv2 / n2
            elif hit_label is not None:
                hx, hy = hit_label
                cx = (bbox.xmin + bbox.xmax) / 2.0
                cy = (bbox.ymin + bbox.ymax) / 2.0
                dxv2, dyv2 = cx - hx, cy - hy
                n2 = np.hypot(dxv2, dyv2)
                if n2 > 1e-6:
                    ux, uy = dxv2 / n2, dyv2 / n2

            lx_disp += ux * STEP_PX
            ly_disp += uy * STEP_PX
            new_xy = inv_trans.transform((lx_disp, ly_disp))
            txt.set_position(tuple(new_xy))
            moved_any = True

        # refresh this label's bbox snapshot for subsequent labels in same pass
        label_bboxes[m] = _bbox_disp(txt)

    if not moved_any:
        break

ax.set_axis_off()
ax.set_title(TITLE, loc="center", fontsize=15, fontweight="bold", pad=22)
# Sub-title with reaction / metabolite counts
ax.text(0.5, 1.015,
        f"{len(net)} reactions \u2022 {G.number_of_nodes()} metabolites",
        transform=ax.transAxes, ha='center', va='bottom',
        fontsize=12, fontweight='bold', color='#222222')

# ---------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------
# Node-type legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Seed (nutrient)',
           markerfacecolor=COLOR_SEED, markeredgecolor='black', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Intermediate metabolite',
           markerfacecolor=COLOR_INTER, markeredgecolor='black', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Core metabolite',
           markerfacecolor=COLOR_CORE, markeredgecolor='black', markersize=10),
]
leg1 = ax.legend(handles=legend_elements, loc='upper left',
                 bbox_to_anchor=(-0.02, 0.47), frameon=False,
                 title="Legend", title_fontsize=12, fontsize=11,
                 prop={'weight': 'bold', 'size': 11})
leg1.get_title().set_fontweight('bold')
ax.add_artist(leg1)

# Degree-size legend
deg_values = [1, 2, 3, 4, 5, 6]
deg_handles = [
    Line2D([0], [0], marker='o', color='w',
           markerfacecolor='black', markeredgecolor='black',
           markersize=np.sqrt(SIZE_BASE + SIZE_STEP * (d - 1)) * 0.55,
           label=str(d))
    for d in deg_values
]
leg2 = ax.legend(handles=deg_handles, loc='upper left',
                 bbox_to_anchor=(-0.02, 1.0), frameon=False,
                 title="Degree", title_fontsize=12, fontsize=11,
                 labelspacing=0.9, borderpad=0.6,
                 prop={'weight': 'bold', 'size': 11})
leg2.get_title().set_fontweight('bold')
ax.add_artist(leg2)

# Growth-medium legend (currency + nutrients listed by name)
medium_kegg = list(Currency) + list(nutrientSet)
medium_names = []
seen = set()
for m in medium_kegg:
    name = met_name(m)
    if name in seen:
        continue
    seen.add(name)
    medium_names.append(name)

medium_text = "Growth medium\n" + "\n".join(medium_names)
ax.text(0.985, 0.985, medium_text,
        transform=ax.transAxes, fontsize=10.5, fontweight='bold',
        va='top', ha='right', multialignment='left',
        bbox=dict(boxstyle="round,pad=0.4", fc="white",
                  ec="black", lw=0.8))

plt.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.04)

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
print(f"Saved figure to {OUT_PATH}")
plt.show()
