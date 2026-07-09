"""
Sweep SAH BVH (2-way) for 128 circles (radius=1), then visualise:
  Left  – spatial view: circles + BVH bounding boxes coloured by depth
  Right – binary tree diagram
"""

import matplotlib
matplotlib.use('Agg')   # non-interactive, no window needed
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from dataclasses import dataclass, field
from typing import List

# ─── AABB (2-D) ──────────────────────────────────────────────────────────────

class AABB:
    def __init__(self):
        self.lo = np.array([ 1e30,  1e30])
        self.hi = np.array([-1e30, -1e30])

    def expand_box(self, b):
        self.lo = np.minimum(self.lo, b.lo)
        self.hi = np.maximum(self.hi, b.hi)

    def half_perim(self):           # 2-D analogue of half surface area
        d = self.hi - self.lo
        return d[0] + d[1]

    def centroid(self):
        return (self.lo + self.hi) * 0.5

    def copy(self):
        b = AABB()
        b.lo, b.hi = self.lo.copy(), self.hi.copy()
        return b

    @staticmethod
    def from_circle(cx, cy, r):
        b = AABB()
        b.lo = np.array([cx - r, cy - r])
        b.hi = np.array([cx + r, cy + r])
        return b

# ─── BVH node ─────────────────────────────────────────────────────────────────

class BVHNode:
    def __init__(self, depth=0):
        self.box: AABB = None
        self.left  = -1
        self.right = -1
        self.prim_indices: List[int] = []
        self.depth = depth

    def is_leaf(self):
        return len(self.prim_indices) > 0

# ─── Sweep SAH builder ────────────────────────────────────────────────────────

class BVHBuilder:
    def __init__(self, c_trav=1.0, c_isect=1.0, max_leaf=4):
        self.c_trav  = c_trav
        self.c_isect = c_isect
        self.max_leaf = max_leaf
        self.nodes: List[BVHNode] = []

    def build(self, boxes: List[AABB]):
        prims = [(b.copy(), i) for i, b in enumerate(boxes)]
        self.nodes = []
        self._rec(prims, 0)

    def _rec(self, prims, depth):
        idx  = len(self.nodes)
        node = BVHNode(depth=depth)
        self.nodes.append(node)

        # union box
        box = AABB()
        for b, _ in prims:
            box.expand_box(b)
        node.box = box
        n = len(prims)

        # ── leaf? ──────────────────────────────────────────────────────────
        if n <= self.max_leaf:
            node.prim_indices = [i for _, i in prims]
            return idx

        # ── find best SAH split ───────────────────────────────────────────
        best_axis  = -1
        best_split = -1
        best_cost  = self.c_isect * n * box.half_perim()   # no-split cost
        best_order = None

        for axis in range(2):
            sp = sorted(prims, key=lambda p: p[0].centroid()[axis])

            # prefix boxes
            L = []
            b = AABB()
            for p, _ in sp:
                b.expand_box(p); L.append(b.copy())

            # suffix boxes
            R = [None] * n
            b = AABB()
            for i in range(n-1, -1, -1):
                b.expand_box(sp[i][0]); R[i] = b.copy()

            inv = 1.0 / box.half_perim() if box.half_perim() > 1e-9 else 0.0

            for s in range(n - 1):
                nl, nr = s + 1, n - s - 1
                cost = (self.c_trav +
                        self.c_isect * inv *
                        (L[s].half_perim() * nl + R[s+1].half_perim() * nr))
                if cost < best_cost:
                    best_cost, best_axis, best_split, best_order = cost, axis, s, sp

        if best_axis < 0:
            node.prim_indices = [i for _, i in prims]
            return idx

        mid = best_split + 1
        l = self._rec(best_order[:mid], depth + 1)
        r = self._rec(best_order[mid:], depth + 1)
        self.nodes[idx].left  = l
        self.nodes[idx].right = r
        return idx

# ─── Build ────────────────────────────────────────────────────────────────────

rng = np.random.default_rng(42)
N, R = 128, 1.0
centers = rng.uniform(-12, 12, (N, 2))

boxes = [AABB.from_circle(cx, cy, R) for cx, cy in centers]

builder = BVHBuilder(max_leaf=4)
builder.build(boxes)
nodes = builder.nodes

max_depth = max(nd.depth for nd in nodes)
n_leaf    = sum(1 for nd in nodes if nd.is_leaf())
n_inner   = len(nodes) - n_leaf
print(f"Nodes total={len(nodes)}  inner={n_inner}  leaf={n_leaf}  max_depth={max_depth}")

# ─── Tree layout (x = leaf order, y = -depth) ────────────────────────────────

nx = [0.0] * len(nodes)
ny = [0.0] * len(nodes)
leaf_ctr = [0]

def layout(i, d):
    nd = nodes[i]
    ny[i] = -d
    if nd.is_leaf():
        nx[i] = float(leaf_ctr[0]); leaf_ctr[0] += 1
    else:
        layout(nd.left,  d + 1)
        layout(nd.right, d + 1)
        nx[i] = (nx[nd.left] + nx[nd.right]) / 2.0

layout(0, 0)

# ─── Plot ─────────────────────────────────────────────────────────────────────

cmap  = plt.cm.plasma
norm  = plt.Normalize(0, max_depth)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10))
fig.patch.set_facecolor('#1a1a2e')
for ax in (ax1, ax2):
    ax.set_facecolor('#16213e')
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('#444')

# ── Spatial view ──────────────────────────────────────────────────────────────
ax1.set_aspect('equal')
ax1.set_title('BVH Spatial View  (128 circles, r=1)', color='white', fontsize=13, pad=10)

# BVH boxes, back to front (deepest depth first so root is on top)
for nd in sorted(nodes, key=lambda x: -x.depth):
    if nd.is_leaf():
        continue
    color = cmap(norm(nd.depth))
    lw    = max(0.3, 1.8 - nd.depth * 0.22)
    b     = nd.box
    rect  = Rectangle((b.lo[0], b.lo[1]), b.hi[0]-b.lo[0], b.hi[1]-b.lo[1],
                       fill=False, edgecolor=color, linewidth=lw, alpha=0.75, zorder=2)
    ax1.add_patch(rect)

# Circles
for cx, cy in centers:
    c = Circle((cx, cy), R, facecolor='#4fc3f7', edgecolor='#0288d1',
                linewidth=0.6, alpha=0.55, zorder=3)
    ax1.add_patch(c)
    ax1.plot(cx, cy, '.', color='#e0e0e0', markersize=1.5, zorder=4)

ax1.set_xlim(centers[:,0].min()-R-1, centers[:,0].max()+R+1)
ax1.set_ylim(centers[:,1].min()-R-1, centers[:,1].max()+R+1)
ax1.set_xlabel('X', color='white'); ax1.set_ylabel('Y', color='white')

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cb = plt.colorbar(sm, ax=ax1, fraction=0.03, pad=0.02)
cb.set_label('BVH depth', color='white')
cb.ax.yaxis.set_tick_params(color='white')
plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')

# ── Tree diagram ──────────────────────────────────────────────────────────────
ax2.set_title(f'BVH Binary Tree  ({len(nodes)} nodes, max depth={max_depth})',
              color='white', fontsize=13, pad=10)
ax2.axis('off')

# edges
for i, nd in enumerate(nodes):
    if not nd.is_leaf():
        for ch in (nd.left, nd.right):
            ax2.plot([nx[i], nx[ch]], [ny[i], ny[ch]],
                     color='#555577', linewidth=0.6, zorder=1)

# nodes
NODE_R = 0.28
for i, nd in enumerate(nodes):
    color = cmap(norm(nd.depth))
    if nd.is_leaf():
        s = NODE_R * 1.6
        rect = FancyBboxPatch((nx[i]-s/2, ny[i]-s/2), s, s,
                               boxstyle="round,pad=0.05",
                               facecolor=color, edgecolor='white',
                               linewidth=0.4, zorder=3)
        ax2.add_patch(rect)
        ax2.text(nx[i], ny[i], str(len(nd.prim_indices)),
                 ha='center', va='center', fontsize=4.5,
                 color='white', fontweight='bold', zorder=4)
    else:
        circ = Circle((nx[i], ny[i]), NODE_R,
                      facecolor=color, edgecolor='white',
                      linewidth=0.4, zorder=3)
        ax2.add_patch(circ)

# depth labels on left
for d in range(max_depth + 1):
    ax2.text(-1.5, -d, f'd={d}', ha='right', va='center',
             fontsize=7, color='#aaa')

# legend patches
legend_elems = [
    mpatches.Patch(color=cmap(norm(d)), label=f'depth {d}')
    for d in range(min(max_depth+1, 6))
]
ax2.legend(handles=legend_elems, loc='lower right',
           fontsize=7, facecolor='#1a1a2e', labelcolor='white',
           edgecolor='#555', framealpha=0.8)

ax2.set_xlim(-3, leaf_ctr[0] + 1)
ax2.set_ylim(-max_depth - 1, 1)

plt.tight_layout()
out = r'C:\Users\ppzh2\OneDrive\Documents\groupmeeting\bvh_result.png'
plt.savefig(out, dpi=160, bbox_inches='tight', facecolor=fig.get_facecolor())
print(f"Saved: {out}")
