"""
4-way BVH comparison — C++ computes, Python visualises.
Three test cases: N=128 / 1024 / 102400
"""

import json, os, subprocess, sys, csv
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Polygon as MplPolygon
import matplotlib.font_manager as fm

BASE_DIR = r'C:\Users\ppzh2\OneDrive\Documents\groupmeeting'
CPP_SRC  = os.path.join(BASE_DIR, 'bvh_compute.cpp')
CPP_EXE  = os.path.join(BASE_DIR, 'bvh_compute.exe')
VCVARS   = (r'C:\Program Files\Microsoft Visual Studio\2022\Community'
            r'\VC\Auxiliary\Build\vcvars64.bat')

# threshold above which we skip circle drawing and full tree diagram
LARGE_N = 1024

for _f in ['Meiryo','Yu Gothic','MS Gothic','BIZ UDGothic','SimHei']:
    if any(ff.name==_f for ff in fm.fontManager.ttflist):
        plt.rcParams['font.family']=_f; break

# ─── AABB / Node ──────────────────────────────────────────────────────────────
class AABB:
    def __init__(self,lo=None,hi=None):
        self.lo = lo if lo is not None else np.array([ 1e30, 1e30])
        self.hi = hi if hi is not None else np.array([-1e30,-1e30])

class Node:
    def __init__(self):
        self.box=AABB(); self.children=[]; self.prims=[]; self.depth=0
    def is_leaf(self): return len(self.prims)>0

def json_to_tree(jnodes):
    nodes=[]
    for jn in jnodes:
        n=Node(); b=jn['box']
        n.box=AABB(np.array([b[0],b[1]]),np.array([b[2],b[3]]))
        n.prims=jn['prims']; n.depth=jn['depth']; nodes.append(n)
    for i,jn in enumerate(jnodes):
        for ci in jn['children']: nodes[i].children.append(nodes[ci])
    return nodes[0]

# ─── Compile & run ────────────────────────────────────────────────────────────
def compile_cpp():
    need=(not os.path.exists(CPP_EXE) or
          os.path.getmtime(CPP_SRC)>os.path.getmtime(CPP_EXE))
    if not need: print("C++ binary up-to-date."); return
    print("Compiling bvh_compute.cpp ...")
    cmd=(f'call "{VCVARS}" >nul 2>&1 && '
         f'cl /O2 /std:c++17 /EHsc /nologo "{CPP_SRC}" /Fe:"{CPP_EXE}"')
    r=subprocess.run(cmd,shell=True,capture_output=True,text=True,cwd=BASE_DIR)
    print(r.stdout.strip())
    if r.returncode!=0: print("Error:\n",r.stderr); sys.exit(1)
    print("Compiled OK ->",CPP_EXE)

def run_cpp():
    r=subprocess.run([CPP_EXE],capture_output=True,text=True,cwd=BASE_DIR)
    print(r.stderr.strip())
    if r.returncode!=0: print("Runtime error:",r.stderr); sys.exit(1)
    return json.loads(r.stdout)

# ─── Tree layout ──────────────────────────────────────────────────────────────
def layout_tree(root, max_draw_depth=None):
    """Returns pos dict {id(node): (x, y)}.
    If max_draw_depth set, only layout nodes up to that depth."""
    pos={}; ctr=[0]
    def _lay(nd,d):
        if max_draw_depth is not None and d>max_draw_depth: return
        if nd.is_leaf() or (max_draw_depth is not None and d==max_draw_depth):
            pos[id(nd)]=(float(ctr[0]),-d); ctr[0]+=1
        else:
            for ch in nd.children: _lay(ch,d+1)
            vis=[id(ch) for ch in nd.children if id(ch) in pos]
            if vis: pos[id(nd)]=(sum(pos[v][0] for v in vis)/len(vis),-d)
            else:   pos[id(nd)]=(float(ctr[0]),-d); ctr[0]+=1
    _lay(root,0); return pos

# ─── Draw tree panel ──────────────────────────────────────────────────────────
def draw_tree(ax, root, title, dt_ms, sah, n_nodes, n_leaves, max_d,
              cmap, max_depth_global, max_draw_depth=None):
    pos  = layout_tree(root, max_draw_depth)
    norm = plt.Normalize(0, max_depth_global)
    ax.set_facecolor('#16213e')
    ax.set_title(title, color='white', fontsize=10, pad=6)
    ax.axis('off')

    # edges
    stack=[root]
    while stack:
        nd=stack.pop()
        if id(nd) not in pos: continue
        x0,y0=pos[id(nd)]
        for ch in nd.children:
            if id(ch) in pos:
                x1,y1=pos[id(ch)]
                ax.plot([x0,x1],[y0,y1],color='#556',linewidth=0.5,zorder=1)
                stack.append(ch)

    # nodes
    RN=0.25; stack=[root]
    while stack:
        nd=stack.pop()
        if id(nd) not in pos: continue
        x,y=pos[id(nd)]; color=cmap(norm(nd.depth))
        is_cut = (max_draw_depth is not None
                  and nd.depth==max_draw_depth and not nd.is_leaf())
        if nd.is_leaf() or is_cut:
            s=RN*1.8
            ax.add_patch(FancyBboxPatch((x-s/2,y-s/2),s,s,
                boxstyle="round,pad=0.04",facecolor=color,edgecolor='white',
                linewidth=0.35,zorder=3))
            label=f"{len(nd.prims)}" if nd.is_leaf() else "…"
            ax.text(x,y,label,ha='center',va='center',
                    fontsize=4,color='white',fontweight='bold',zorder=4)
        else:
            ax.add_patch(Circle((x,y),RN,facecolor=color,edgecolor='white',
                linewidth=0.35,zorder=3))
            for ch in nd.children: stack.append(ch)

    cut_note = f"  (top {max_draw_depth} levels shown)" if max_draw_depth else ""
    info=(f"Build : {dt_ms:.4f} ms\n"
          f"SAH   : {sah:.4f}\n"
          f"Nodes : {n_nodes}  Leaves: {n_leaves}\n"
          f"Depth : {max_d}{cut_note}")
    ax.text(0.02,0.02,info,transform=ax.transAxes,fontsize=7,color='#ccc',
        va='bottom',fontfamily='monospace',
        bbox=dict(boxstyle='round',facecolor='#0d0d1e',alpha=0.7))
    xs=[p[0] for p in pos.values()]; ys=[p[1] for p in pos.values()]
    if xs: ax.set_xlim(min(xs)-1,max(xs)+1); ax.set_ylim(min(ys)-1,max(ys)+1)

# ─── Draw spatial panel ───────────────────────────────────────────────────────
def draw_spatial(ax, root, prim_data, prim_type, R_val, N, title, cmap, max_depth_global,
                 max_box_depth=None):
    ax.set_facecolor('#16213e'); ax.set_aspect('equal')
    ax.set_title(title,color='white',fontsize=10,pad=6)
    ax.tick_params(colors='white')
    for sp in ax.spines.values(): sp.set_edgecolor('#444')
    norm=plt.Normalize(0,max_depth_global)

    all_nodes=[]
    stack=[root]
    while stack:
        nd=stack.pop(); all_nodes.append(nd)
        if max_box_depth is None or nd.depth < max_box_depth:
            for ch in nd.children: stack.append(ch)

    for nd in sorted(all_nodes, key=lambda x: -x.depth):
        color=cmap(norm(nd.depth))
        if nd.is_leaf():
            lw,ls,alpha,z=0.9,'--',0.9,4
        else:
            lw,ls,alpha,z=max(0.4,1.8-nd.depth*0.22),'-',0.75,2
        b=nd.box
        ax.add_patch(Rectangle((b.lo[0],b.lo[1]),b.hi[0]-b.lo[0],b.hi[1]-b.lo[1],
            fill=False,edgecolor=color,linewidth=lw,linestyle=ls,alpha=alpha,zorder=z))

    # Draw primitives
    if prim_type == 'triangle':
        tri_arr = np.array(prim_data, dtype=np.float32)  # (N, 3, 2)
        if N <= LARGE_N:
            for verts in tri_arr:
                ax.add_patch(MplPolygon(verts, closed=True,
                    facecolor='#4fc3f7', edgecolor='#0288d1',
                    linewidth=0.3, alpha=0.45, zorder=3))
        else:
            cents = tri_arr.mean(axis=1)   # (N, 2) centroids
            ax.scatter(cents[:,0], cents[:,1], s=0.2, c='#4fc3f7', alpha=0.3, zorder=3)
        all_v = tri_arr.reshape(-1, 2)
        ax.set_xlim(all_v[:,0].min()-1, all_v[:,0].max()+1)
        ax.set_ylim(all_v[:,1].min()-1, all_v[:,1].max()+1)
        if N > LARGE_N:
            ax.text(0.02,0.98,f"N={N:,}  (dots=triangle centroids)",
                    transform=ax.transAxes,color='#aaa',fontsize=7,va='top')
    else:
        centers = np.array(prim_data, dtype=np.float32)
        if N <= LARGE_N:
            for cx,cy in centers:
                ax.add_patch(Circle((cx,cy),R_val,facecolor='#4fc3f7',
                    edgecolor='#0288d1',linewidth=0.3,alpha=0.45,zorder=3))
        else:
            ax.scatter(centers[:,0],centers[:,1],s=0.2,c='#4fc3f7',alpha=0.3,zorder=3)
        xpad=R_val+1
        ax.set_xlim(centers[:,0].min()-xpad, centers[:,0].max()+xpad)
        ax.set_ylim(centers[:,1].min()-xpad, centers[:,1].max()+xpad)
        if N > LARGE_N:
            ax.text(0.02,0.98,f"N={N:,}  (dots=circle centres)",
                    transform=ax.transAxes,color='#aaa',fontsize=7,va='top')
    ax.set_xlabel('X',color='white',fontsize=8)
    ax.set_ylabel('Y',color='white',fontsize=8)

# ─── Slug map ─────────────────────────────────────────────────────────────────
SLUG = {
    'Median Split':         '0_median',
    '2-way SAH Sweep':      '1_2way',
    'Binned SAH (B=16)':    '1b_binned_16',
    'Binned Collapse k=2':  '1c_binned_collapse',
    'Collapse k=2':         '2_collapse_k2',
    'A Independent':        '3_A_independent',
    'B Hierarchical':       '4_B_hierarchical',
    'C Exhaustive 2D':      '5_C_exhaustive',
    'D All-Axis':           '6_D_all_axis',
    'Binned+A (T=32)':          '7_binned_A_32',
    'Binned+SAH (T=32)':        '8_binned_SAH_32',
    'Binned4+A (T=32)':         '9_binned4_A_32',
    'Binned4+Coll (T=32)':      '9c_binned4_coll_32',
    'Binned+SAH+Coll (T=32)':   '9b_binned_SAH_coll_32',
}

# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════
compile_cpp()
data = run_cpp()
BG   = '#1a1a2e'
cmap = plt.cm.plasma

# ─── Timestamped output folder ────────────────────────────────────────────────
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
OUT_DIR   = os.path.join(BASE_DIR, f'run_{timestamp}')
os.makedirs(OUT_DIR, exist_ok=True)
print(f"\nOutput folder: {OUT_DIR}\n")

all_results = []   # for summary table: list of (N, strategy_name, dt_ms, sah, ...)

for tc in data['test_cases']:
    N         = tc['N']
    R_val     = float(tc['R'])
    prim_type = tc.get('type', 'circle')
    prim_data = tc['triangles'] if prim_type == 'triangle' else tc['circles']
    type_tag  = '_tri' if prim_type == 'triangle' else ''
    strats    = tc['strategies']

    max_depth_global = max(st['max_depth'] for st in strats)

    # How many tree levels to draw
    max_draw_depth = None if N <= LARGE_N else 5
    max_box_depth  = None if N <= LARGE_N else 6

    print(f"\n--- Generating figures for N={N} ({prim_type}) ---")

    for st in strats:
        name  = st['name']
        slug  = SLUG.get(name, name.replace(' ','_'))
        root  = json_to_tree(st['nodes'])
        dt_ms = st['build_time_ms']
        sah   = st['sah_cost']
        nd_c  = st['node_count']
        lf_c  = st['leaf_count']
        md    = st['max_depth']

        all_results.append((N, prim_type, name, dt_ms, sah, nd_c, lf_c, md))

        fig, (ax_sp, ax_tr) = plt.subplots(1, 2, figsize=(16, 8))
        fig.patch.set_facecolor(BG)
        prim_label = 'triangles' if prim_type == 'triangle' else 'circles'
        fig.text(0.5,0.98,f"BVH  [N={N:,}  {prim_label}]  --  {name}",
            ha='center',va='top',color='white',fontsize=13,fontweight='bold')
        fig.text(0.5,0.95,
            "Left: spatial AABB overlay (solid=inner, dashed=leaf, colour=depth).  "
            "Right: tree structure.  C++ computation.",
            ha='center',va='top',color='#aaa',fontsize=9)

        draw_spatial(ax_sp, root, prim_data, prim_type, R_val, N,
                     f"Spatial View  (N={N:,}  {prim_label})", cmap, max_depth_global, max_box_depth)
        draw_tree(ax_tr, root, "Tree Structure",
                  dt_ms, sah, nd_c, lf_c, md,
                  cmap, max_depth_global, max_draw_depth)

        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,max_depth_global))
        cb = fig.colorbar(sm,ax=[ax_sp,ax_tr],fraction=0.018,pad=0.01,aspect=30)
        cb.set_label('BVH depth',color='white',fontsize=9)
        cb.ax.yaxis.set_tick_params(color='white')
        plt.setp(cb.ax.yaxis.get_ticklabels(),color='white')

        plt.subplots_adjust(top=0.93,bottom=0.04,left=0.04,right=0.92,wspace=0.1)
        out = os.path.join(OUT_DIR, f'bvh4_N{N}{type_tag}_{slug}.png')
        plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
        plt.close()
        print(f"  Saved -> {out}")

    # ── Per-test-case table ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 2.5 + 0.5*len(strats)))
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG); ax.axis('off')
    fig.text(0.5,0.97,f"BVH Strategy Comparison  (N={N:,}, {prim_type}, r={R_val}, range={tc['range']})",
        ha='center',va='top',color='white',fontsize=12,fontweight='bold')

    col_labels=['Strategy','Build Time','SAH Cost','Total Nodes','Leaf Nodes','Max Depth']
    rows=[[st['name'], f"{st['build_time_ms']:.4f} ms", f"{st['sah_cost']:.4f}",
           str(st['node_count']),str(st['leaf_count']),str(st['max_depth'])]
          for st in strats]

    tbl=ax.table(cellText=rows,colLabels=col_labels,loc='center',cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1,2.0)

    nr=len(rows)
    best={1:min(range(nr),key=lambda i:strats[i]['build_time_ms']),
          2:min(range(nr),key=lambda i:strats[i]['sah_cost']),
          3:min(range(nr),key=lambda i:strats[i]['node_count']),
          5:min(range(nr),key=lambda i:strats[i]['max_depth'])}
    for j in range(len(col_labels)):
        tbl[0,j].set_facecolor('#2e2e5e')
        tbl[0,j].set_text_props(color='white',fontweight='bold')
    for i in range(nr):
        for j in range(len(col_labels)):
            cell=tbl[i+1,j]
            cell.set_facecolor('#1a1a3e' if i%2==0 else '#16213e')
            cell.set_text_props(color='white')
            if j in best and best[j]==i: cell.set_facecolor('#1a4a2e')

    out=os.path.join(OUT_DIR,f'bvh4_N{N}{type_tag}_table.png')
    plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
    plt.close(); print(f"  Saved -> {out}")

    # ── Per-test-case CSV ──────────────────────────────────────────────────
    csv_path = os.path.join(OUT_DIR, f'bvh4_N{N}{type_tag}_table.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['Strategy','Build Time (ms)','SAH Cost',
                    'Total Nodes','Leaf Nodes','Max Depth'])
        for st in strats:
            w.writerow([st['name'],
                        f"{st['build_time_ms']:.4f}",
                        f"{st['sah_cost']:.4f}",
                        st['node_count'], st['leaf_count'], st['max_depth']])
    print(f"  Saved -> {csv_path}")

# ─── Summary table: build time scaling across N ───────────────────────────────
strat_names = list(dict.fromkeys(r[2] for r in all_results))
NP_values   = list(dict.fromkeys((r[0],r[1]) for r in all_results))  # (N, prim_type)

# build time / SAH lookup keyed by (N, prim_type, strategy)
lookup_dt  = {(r[0],r[1],r[2]): r[3] for r in all_results}
lookup_sah = {(r[0],r[1],r[2]): r[4] for r in all_results}

for metric, lookup, label, fmt in [
    ('Build Time (ms)', lookup_dt,  'bvh4_summary_buildtime.png', lambda v: f"{v:.4f}"),
    ('SAH Cost',        lookup_sah, 'bvh4_summary_sah.png',       lambda v: f"{v:.2f}"),
]:
    fig, ax = plt.subplots(figsize=(14, 2 + 0.55*len(NP_values)))
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG); ax.axis('off')
    fig.text(0.5,0.97,f"Scaling Comparison — {metric}",
        ha='center',va='top',color='white',fontsize=12,fontweight='bold')
    fig.text(0.5,0.90,"Green = best in row.  N/A = strategy skipped (too slow).",
        ha='center',va='top',color='#aaa',fontsize=9)

    col_labels=['N (prim)'] + strat_names
    rows=[]
    for N, pt in NP_values:
        row=[f"{N:,} ({pt[:3]})"]
        for sn in strat_names:
            v=lookup.get((N,pt,sn),None)
            row.append(fmt(v) if v is not None else "N/A")
        rows.append(row)

    tbl=ax.table(cellText=rows,colLabels=col_labels,loc='center',cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1,2.0)
    for j in range(len(col_labels)):
        tbl[0,j].set_facecolor('#2e2e5e')
        tbl[0,j].set_text_props(color='white',fontweight='bold')
    for i,(N,pt) in enumerate(NP_values):
        vals=[(j+1, lookup.get((N,pt,sn),None)) for j,sn in enumerate(strat_names)
              if lookup.get((N,pt,sn)) is not None]
        best_j=min(vals,key=lambda x:x[1])[0] if vals else -1
        for j in range(len(col_labels)):
            cell=tbl[i+1,j]
            cell.set_facecolor('#1a1a3e' if i%2==0 else '#16213e')
            cell.set_text_props(color='white')
            if j==best_j: cell.set_facecolor('#1a4a2e')

    out=os.path.join(OUT_DIR,label)
    plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
    plt.close(); print(f"Saved -> {out}")

    # ── Summary CSV ───────────────────────────────────────────────────────
    csv_label = label.replace('.png', '.csv')
    csv_path  = os.path.join(OUT_DIR, csv_label)
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['N (prim)'] + strat_names)
        for N, pt in NP_values:
            row = [f"{N:,} ({pt[:3]})"]
            for sn in strat_names:
                v = lookup.get((N, pt, sn), None)
                row.append(fmt(v) if v is not None else "N/A")
            w.writerow(row)
    print(f"Saved -> {csv_path}")

print("\nAll done.")
