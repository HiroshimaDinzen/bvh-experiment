// bvh_compute.cpp
// Multi-test BVH comparison: N=128/1024/102400, 5 strategies.
// Outputs JSON to stdout, progress to stderr.
// Compile (MSVC): cl /O2 /std:c++17 /EHsc bvh_compute.cpp /Fe:bvh_compute.exe
// Compile (g++):  g++ -O2 -std=c++17 bvh_compute.cpp -o bvh_compute

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>

static const float C_TRAV        = 1.0f;
static const float C_ISECT       = 1.0f;
static const int   MAX_LEAF      = 4;
static const int   HYBRID_THRESH = 32;  // switch from Binned to precise method below this count

// ─── AABB 2-D ─────────────────────────────────────────────────────────────────
struct AABB {
    float lo[2] = {  1e30f,  1e30f };
    float hi[2] = { -1e30f, -1e30f };
    void expand(const AABB& b) {
        lo[0]=std::min(lo[0],b.lo[0]); lo[1]=std::min(lo[1],b.lo[1]);
        hi[0]=std::max(hi[0],b.hi[0]); hi[1]=std::max(hi[1],b.hi[1]);
    }
    float half_perim()     const { return (hi[0]-lo[0])+(hi[1]-lo[1]); }
    float centroid(int ax) const { return (lo[ax]+hi[ax])*0.5f; }
    static AABB circle(float cx,float cy,float r) {
        AABB b; b.lo[0]=cx-r;b.lo[1]=cy-r;b.hi[0]=cx+r;b.hi[1]=cy+r; return b;
    }
};

struct Prim { AABB box; int idx; };
using  Prims = std::vector<Prim>;

AABB union_box(const Prims& ps){
    AABB b; for(auto& p:ps) b.expand(p.box); return b;
}

// ─── BVH Node ─────────────────────────────────────────────────────────────────
struct Node {
    AABB             box;
    std::vector<int> children;
    std::vector<int> prims;
    int              depth=0;
    bool is_leaf() const { return !prims.empty(); }
};

// ─── SAH cost (whole tree) ────────────────────────────────────────────────────
double sah_cost(const std::vector<Node>& nodes){
    float ra=nodes[0].box.half_perim();
    if(ra<1e-9f) return 0.0;
    double c=0.0;
    for(auto& n:nodes){
        float w=n.box.half_perim()/ra;
        c += n.is_leaf() ? C_ISECT*(int)n.prims.size()*w : C_TRAV*w;
    }
    return c;
}

struct Stats { int nodes,leaves,max_depth,prim_count; };
Stats tree_stats(const std::vector<Node>& nodes){
    Stats s{0,0,0,0};
    for(auto& n:nodes){
        ++s.nodes;
        if(n.depth>s.max_depth) s.max_depth=n.depth;
        if(n.is_leaf()){ ++s.leaves; s.prim_count+=(int)n.prims.size(); }
    }
    return s;
}

// ─── 1-D SAH sweep ────────────────────────────────────────────────────────────
struct Split { int s; float cost; Prims sorted; };

Split best_1d(const Prims& ps, int axis){
    Prims sp=ps;
    std::sort(sp.begin(),sp.end(),[axis](const Prim& a,const Prim& b){
        return a.box.centroid(axis)<b.box.centroid(axis);
    });
    int n=(int)sp.size();
    float pa=union_box(sp).half_perim();
    float inv=pa>1e-9f?1.f/pa:0.f;
    std::vector<float> L(n),R(n);
    { AABB b; for(int i=0;i<n;i++)    {b.expand(sp[i].box);L[i]=b.half_perim();} }
    { AABB b; for(int i=n-1;i>=0;i--) {b.expand(sp[i].box);R[i]=b.half_perim();} }
    float best=C_ISECT*n; int bs=-1;
    for(int s=0;s<n-1;s++){
        float c=C_TRAV+C_ISECT*inv*(L[s]*(s+1)+R[s+1]*(n-s-1));
        if(c<best){best=c;bs=s;}
    }
    return {bs,best,std::move(sp)};
}

Split best_any(const Prims& ps){
    Split best{-1,C_ISECT*(float)ps.size(),ps};
    for(int ax=0;ax<2;ax++){
        auto r=best_1d(ps,ax);
        if(r.s>=0&&r.cost<best.cost) best=std::move(r);
    }
    return best;
}

static void  make_leaf(std::vector<Node>& ns,int i,const Prims& ps){for(auto& p:ps)ns[i].prims.push_back(p.idx);}
static int   new_node (std::vector<Node>& ns,const Prims& ps,int d){int i=(int)ns.size();ns.push_back(Node());ns[i].box=union_box(ps);ns[i].depth=d;return i;}
static Prims slice    (const Prims& sp,int a,int b){return Prims(sp.begin()+a,sp.begin()+b);}

// ════════════════════════════════════════════════════════════════════════════
// Build functions
// ════════════════════════════════════════════════════════════════════════════
int build_median(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    float sp[2]={};
    for(int ax=0;ax<2;ax++){
        float lo=1e30f,hi=-1e30f;
        for(auto& p:ps){float c=p.box.centroid(ax);lo=std::min(lo,c);hi=std::max(hi,c);}
        sp[ax]=hi-lo;
    }
    int axis=sp[0]>=sp[1]?0:1;
    Prims sorted=ps;
    std::sort(sorted.begin(),sorted.end(),[axis](const Prim& a,const Prim& b){
        return a.box.centroid(axis)<b.box.centroid(axis);
    });
    int mid=n/2;
    int l=build_median(ns,slice(sorted,0,mid),d+1);
    int r=build_median(ns,slice(sorted,mid,n),d+1);
    ns[idx].children={l,r}; return idx;
}

int build_2way(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    auto r=best_any(ps);
    if(r.s<0){make_leaf(ns,idx,ps);return idx;}
    int l =build_2way(ns,slice(r.sorted,0,r.s+1),d+1);
    int ri=build_2way(ns,slice(r.sorted,r.s+1,(int)r.sorted.size()),d+1);
    ns[idx].children={l,ri}; return idx;
}

// ─── Binned SAH (2-way, B bins per axis) ──────────────────────────────────────
// Instead of evaluating all N-1 centroid-sorted split positions (sweep SAH),
// divide the centroid range into B equal-width bins and evaluate only the
// B-1 bin-boundary split candidates.  O(N + B) per node vs O(N log N).
int build_binned(std::vector<Node>& ns, const Prims& ps, int d, int n_bins) {
    int idx = new_node(ns, ps, d);
    int n   = (int)ps.size();
    if (n <= MAX_LEAF) { make_leaf(ns, idx, ps); return idx; }

    float pa  = ns[idx].box.half_perim();
    float inv = pa > 1e-9f ? 1.f / pa : 0.f;

    float best_cost = C_ISECT * n;   // no-split cost (leaf)
    int   best_axis = -1, best_bin  = -1;

    struct Bin { AABB box; int cnt = 0; };

    for (int axis = 0; axis < 2; axis++) {
        // centroid range along this axis
        float c_lo = 1e30f, c_hi = -1e30f;
        for (auto& p : ps) {
            float c = p.box.centroid(axis);
            c_lo = std::min(c_lo, c);
            c_hi = std::max(c_hi, c);
        }
        if (c_hi - c_lo < 1e-9f) continue;   // all centroids coincide on this axis

        float scale = n_bins / (c_hi - c_lo);
        std::vector<Bin> bins(n_bins);

        // assign each primitive to its bin
        for (auto& p : ps) {
            int b = (int)((p.box.centroid(axis) - c_lo) * scale);
            b = std::min(b, n_bins - 1);      // clamp last prim to final bin
            bins[b].box.expand(p.box);
            bins[b].cnt++;
        }

        // prefix sweep: L[i] = half-perim of union of bins [0..i], Lc[i] = count
        std::vector<float> L(n_bins); std::vector<int> Lc(n_bins, 0);
        { AABB b; int cnt = 0;
          for (int i = 0; i < n_bins; i++) {
              if (bins[i].cnt > 0) b.expand(bins[i].box);
              cnt += bins[i].cnt;
              L[i] = b.half_perim(); Lc[i] = cnt;
          }
        }
        // suffix sweep: R[i] = half-perim of union of bins [i..B-1], Rc[i] = count
        std::vector<float> R(n_bins); std::vector<int> Rc(n_bins, 0);
        { AABB b; int cnt = 0;
          for (int i = n_bins - 1; i >= 0; i--) {
              if (bins[i].cnt > 0) b.expand(bins[i].box);
              cnt += bins[i].cnt;
              R[i] = b.half_perim(); Rc[i] = cnt;
          }
        }

        // evaluate the B-1 bin-boundary split candidates
        for (int s = 0; s < n_bins - 1; s++) {
            int nl = Lc[s], nr = Rc[s + 1];
            if (nl == 0 || nr == 0) continue;
            float cost = C_TRAV + C_ISECT * inv * (L[s] * nl + R[s + 1] * nr);
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_bin  = s;
            }
        }
    }

    if (best_axis < 0) { make_leaf(ns, idx, ps); return idx; }

    // re-partition primitives on the winning (axis, bin-boundary)
    float c_lo = 1e30f, c_hi = -1e30f;
    for (auto& p : ps) {
        float c = p.box.centroid(best_axis);
        c_lo = std::min(c_lo, c);
        c_hi = std::max(c_hi, c);
    }
    float scale = n_bins / (c_hi - c_lo);

    Prims lp, rp;
    for (auto& p : ps) {
        int b = (int)((p.box.centroid(best_axis) - c_lo) * scale);
        b = std::min(b, n_bins - 1);
        (b <= best_bin ? lp : rp).push_back(p);
    }
    if (lp.empty() || rp.empty()) { make_leaf(ns, idx, ps); return idx; }

    int l = build_binned(ns, lp, d + 1, n_bins);
    int r = build_binned(ns, rp, d + 1, n_bins);
    ns[idx].children = {l, r};
    return idx;
}

int build_A(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    auto rx=best_1d(ps,0),ry=best_1d(ps,1);
    if(rx.s<0&&ry.s<0){make_leaf(ns,idx,ps);return idx;}
    Split *pr; int sec;
    if(rx.s>=0&&(ry.s<0||rx.cost<=ry.cost)){pr=&rx;sec=1;}
    else                                    {pr=&ry;sec=0;}
    Prims L=slice(pr->sorted,0,pr->s+1);
    Prims R=slice(pr->sorted,pr->s+1,(int)pr->sorted.size());
    auto rL=best_1d(L,sec),rR=best_1d(R,sec);
    std::vector<Prims> quads;
    if(rL.s>=0){quads.push_back(slice(rL.sorted,0,rL.s+1));
                quads.push_back(slice(rL.sorted,rL.s+1,(int)rL.sorted.size()));}
    else        quads.push_back(L);
    if(rR.s>=0){quads.push_back(slice(rR.sorted,0,rR.s+1));
                quads.push_back(slice(rR.sorted,rR.s+1,(int)rR.sorted.size()));}
    else        quads.push_back(R);
    quads.erase(std::remove_if(quads.begin(),quads.end(),[](const Prims& q){return q.empty();}),quads.end());
    if(quads.size()==1){make_leaf(ns,idx,ps);return idx;}
    for(auto& q:quads) ns[idx].children.push_back(build_A(ns,q,d+1));
    return idx;
}

int build_B(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    auto r1=best_any(ps);
    if(r1.s<0){make_leaf(ns,idx,ps);return idx;}
    Prims L1=slice(r1.sorted,0,r1.s+1);
    Prims R1=slice(r1.sorted,r1.s+1,(int)r1.sorted.size());
    auto r2=best_any(L1),r3=best_any(R1);
    std::vector<Prims> quads;
    if(r2.s>=0){quads.push_back(slice(r2.sorted,0,r2.s+1));
                quads.push_back(slice(r2.sorted,r2.s+1,(int)r2.sorted.size()));}
    else        quads.push_back(L1);
    if(r3.s>=0){quads.push_back(slice(r3.sorted,0,r3.s+1));
                quads.push_back(slice(r3.sorted,r3.s+1,(int)r3.sorted.size()));}
    else        quads.push_back(R1);
    quads.erase(std::remove_if(quads.begin(),quads.end(),[](const Prims& q){return q.empty();}),quads.end());
    if(quads.size()==1){make_leaf(ns,idx,ps);return idx;}
    for(auto& q:quads) ns[idx].children.push_back(build_B(ns,q,d+1));
    return idx;
}

int build_C(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    float pa=ns[idx].box.half_perim();
    float inv=pa>1e-9f?1.f/pa:0.f;
    Prims spx=ps,spy=ps;
    std::sort(spx.begin(),spx.end(),[](const Prim& a,const Prim& b){return a.box.centroid(0)<b.box.centroid(0);});
    std::sort(spy.begin(),spy.end(),[](const Prim& a,const Prim& b){return a.box.centroid(1)<b.box.centroid(1);});
    std::unordered_map<int,int> xrank;
    for(int i=0;i<n;i++) xrank[spx[i].idx]=i;
    std::vector<int> spy_xr(n);
    for(int j=0;j<n;j++) spy_xr[j]=xrank[spy[j].idx];
    float best=C_ISECT*n; int bsx=-1,bsy=-1;
    for(int sx=0;sx<n-1;sx++){
        for(int sy=0;sy<n-1;sy++){
            AABB bLL,bLR,bRL,bRR; int nLL=0,nLR=0,nRL=0,nRR=0;
            for(int j=0;j<n;j++){
                bool left=spy_xr[j]<=sx,bot=j<=sy;
                if     ( left&& bot){bLL.expand(spy[j].box);nLL++;}
                else if( left&&!bot){bLR.expand(spy[j].box);nLR++;}
                else if(!left&& bot){bRL.expand(spy[j].box);nRL++;}
                else                {bRR.expand(spy[j].box);nRR++;}
            }
            float c=C_TRAV;
            if(nLL>0) c+=C_ISECT*inv*bLL.half_perim()*nLL;
            if(nLR>0) c+=C_ISECT*inv*bLR.half_perim()*nLR;
            if(nRL>0) c+=C_ISECT*inv*bRL.half_perim()*nRL;
            if(nRR>0) c+=C_ISECT*inv*bRR.half_perim()*nRR;
            if(c<best){best=c;bsx=sx;bsy=sy;}
        }
    }
    if(bsx<0){make_leaf(ns,idx,ps);return idx;}
    Prims LL,LR,RL,RR;
    for(int j=0;j<n;j++){
        bool left=spy_xr[j]<=bsx,bot=j<=bsy;
        if     ( left&& bot) LL.push_back(spy[j]);
        else if( left&&!bot) LR.push_back(spy[j]);
        else if(!left&& bot) RL.push_back(spy[j]);
        else                 RR.push_back(spy[j]);
    }
    std::vector<Prims> quads;
    for(auto& q:{LL,LR,RL,RR}) if(!q.empty()) quads.push_back(q);
    if(quads.size()==1){make_leaf(ns,idx,ps);return idx;}
    for(auto& q:quads) ns[idx].children.push_back(build_C(ns,q,d+1));
    return idx;
}

// ─── D All-Axis Exhaustive ────────────────────────────────────────────────────
// Extends strategy C by also enumerating same-axis two-cut candidates:
//   Cross-axis  (sx, sy) → up to 4 children on a 2-D grid  [O(N³)]
//   Same-axis X (s1, s2) → 3 children: Left/Middle/Right along X  [O(N²)]
//   Same-axis Y (s1, s2) → 3 children: Left/Middle/Right along Y  [O(N²)]
// Total per-node: O(N³), heavy — skip for large N.
int build_D(std::vector<Node>& ns, const Prims& ps, int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    float pa=ns[idx].box.half_perim();
    float inv=pa>1e-9f?1.f/pa:0.f;

    // sort by each axis; build X-rank lookup
    Prims spx=ps,spy=ps;
    std::sort(spx.begin(),spx.end(),[](const Prim&a,const Prim&b){return a.box.centroid(0)<b.box.centroid(0);});
    std::sort(spy.begin(),spy.end(),[](const Prim&a,const Prim&b){return a.box.centroid(1)<b.box.centroid(1);});
    std::unordered_map<int,int> xrank;
    for(int i=0;i<n;i++) xrank[spx[i].idx]=i;
    std::vector<int> spy_xr(n);
    for(int j=0;j<n;j++) spy_xr[j]=xrank[spy[j].idx];

    // prefix/suffix AABB objects for same-axis middle computation
    std::vector<AABB> Lx(n),Rx(n),Ly(n),Ry(n);
    { AABB b; for(int i=0;i<n;i++){b.expand(spx[i].box);Lx[i]=b;} }
    { AABB b; for(int i=n-1;i>=0;i--){b.expand(spx[i].box);Rx[i]=b;} }
    { AABB b; for(int i=0;i<n;i++){b.expand(spy[i].box);Ly[i]=b;} }
    { AABB b; for(int i=n-1;i>=0;i--){b.expand(spy[i].box);Ry[i]=b;} }

    float best=C_ISECT*n;
    int btype=-1;             // 0=cross, 1=same-X, 2=same-Y
    int bsx=-1,bsy=-1;        // cross-axis winners
    int bs1=-1,bs2=-1;        // same-axis winners

    // ── cross-axis (sx,sy) pairs  [same as strategy C] ───────────────────────
    for(int sx=0;sx<n-1;sx++){
        for(int sy=0;sy<n-1;sy++){
            AABB qLL,qLR,qRL,qRR; int nLL=0,nLR=0,nRL=0,nRR=0;
            for(int j=0;j<n;j++){
                bool lft=spy_xr[j]<=sx, bot=j<=sy;
                if     (lft&& bot){qLL.expand(spy[j].box);nLL++;}
                else if(lft&&!bot){qLR.expand(spy[j].box);nLR++;}
                else if(!lft&&bot){qRL.expand(spy[j].box);nRL++;}
                else              {qRR.expand(spy[j].box);nRR++;}
            }
            float c=C_TRAV;
            if(nLL>0)c+=C_ISECT*inv*qLL.half_perim()*nLL;
            if(nLR>0)c+=C_ISECT*inv*qLR.half_perim()*nLR;
            if(nRL>0)c+=C_ISECT*inv*qRL.half_perim()*nRL;
            if(nRR>0)c+=C_ISECT*inv*qRR.half_perim()*nRR;
            if(c<best){best=c;btype=0;bsx=sx;bsy=sy;}
        }
    }

    // ── same-axis X: two cuts s1 < s2  → 3 children ─────────────────────────
    // Left  = spx[0..s1]      AABB = Lx[s1]
    // Middle= spx[s1+1..s2]   AABB = incremental (reset each outer iter)
    // Right = spx[s2+1..n-1]  AABB = Rx[s2+1]
    for(int s1=0;s1<n-2;s1++){
        AABB mid;
        for(int s2=s1+1;s2<n-1;s2++){
            mid.expand(spx[s2].box);            // mid covers spx[s1+1..s2]
            int nL=s1+1, nM=s2-s1, nR=n-s2-1;
            float c=C_TRAV;
            c+=C_ISECT*inv*Lx[s1].half_perim()*nL;
            c+=C_ISECT*inv*mid.half_perim()*nM;
            c+=C_ISECT*inv*Rx[s2+1].half_perim()*nR;
            if(c<best){best=c;btype=1;bs1=s1;bs2=s2;}
        }
    }

    // ── same-axis Y: two cuts s1 < s2  → 3 children ─────────────────────────
    for(int s1=0;s1<n-2;s1++){
        AABB mid;
        for(int s2=s1+1;s2<n-1;s2++){
            mid.expand(spy[s2].box);            // mid covers spy[s1+1..s2]
            int nL=s1+1, nM=s2-s1, nR=n-s2-1;
            float c=C_TRAV;
            c+=C_ISECT*inv*Ly[s1].half_perim()*nL;
            c+=C_ISECT*inv*mid.half_perim()*nM;
            c+=C_ISECT*inv*Ry[s2+1].half_perim()*nR;
            if(c<best){best=c;btype=2;bs1=s1;bs2=s2;}
        }
    }

    if(btype<0){make_leaf(ns,idx,ps);return idx;}

    // ── partition ─────────────────────────────────────────────────────────────
    std::vector<Prims> groups;
    if(btype==0){
        // 4-way grid split
        Prims gLL,gLR,gRL,gRR;
        for(int j=0;j<n;j++){
            bool lft=spy_xr[j]<=bsx, bot=j<=bsy;
            if     (lft&& bot) gLL.push_back(spy[j]);
            else if(lft&&!bot) gLR.push_back(spy[j]);
            else if(!lft&&bot) gRL.push_back(spy[j]);
            else               gRR.push_back(spy[j]);
        }
        for(auto& g:{gLL,gLR,gRL,gRR}) if(!g.empty()) groups.push_back(g);
    } else {
        // 3-way same-axis split
        const Prims& sp=(btype==1)?spx:spy;
        Prims gL=slice(sp,0,bs1+1);
        Prims gM=slice(sp,bs1+1,bs2+1);
        Prims gR=slice(sp,bs2+1,n);
        for(auto& g:{gL,gM,gR}) if(!g.empty()) groups.push_back(g);
    }
    if(groups.size()<=1){make_leaf(ns,idx,ps);return idx;}
    for(auto& g:groups) ns[idx].children.push_back(build_D(ns,g,d+1));
    return idx;
}

// ─── Hybrid helpers ───────────────────────────────────────────────────────────
// binned_split_2way: perform one Binned SAH (B=n_bins) 2-way split.
//   pa     = half-perimeter of parent AABB (pre-computed)
//   lp/rp  = output left / right primitive sets
//   returns true iff a valid split was found (both sides non-empty)
static bool binned_split_2way(const Prims& ps, float pa, int n_bins,
                               Prims& lp, Prims& rp)
{
    int n = (int)ps.size();
    float inv = pa > 1e-9f ? 1.f / pa : 0.f;
    float best_cost = C_ISECT * n;
    int   best_axis = -1, best_bin  = -1;
    float best_clo  = 0.f, best_scale = 0.f;

    struct Bin { AABB box; int cnt = 0; };
    for (int axis = 0; axis < 2; axis++) {
        float c_lo = 1e30f, c_hi = -1e30f;
        for (auto& p : ps) { float c = p.box.centroid(axis); c_lo=std::min(c_lo,c); c_hi=std::max(c_hi,c); }
        if (c_hi - c_lo < 1e-9f) continue;
        float scale = n_bins / (c_hi - c_lo);
        std::vector<Bin> bins(n_bins);
        for (auto& p : ps) {
            int b = std::min((int)((p.box.centroid(axis) - c_lo) * scale), n_bins - 1);
            bins[b].box.expand(p.box); bins[b].cnt++;
        }
        std::vector<float> L(n_bins); std::vector<int> Lc(n_bins, 0);
        { AABB b; int c=0; for(int i=0;i<n_bins;i++){ if(bins[i].cnt>0)b.expand(bins[i].box); c+=bins[i].cnt; L[i]=b.half_perim(); Lc[i]=c; } }
        std::vector<float> R(n_bins); std::vector<int> Rc(n_bins, 0);
        { AABB b; int c=0; for(int i=n_bins-1;i>=0;i--){ if(bins[i].cnt>0)b.expand(bins[i].box); c+=bins[i].cnt; R[i]=b.half_perim(); Rc[i]=c; } }
        for (int s = 0; s < n_bins - 1; s++) {
            int nl=Lc[s], nr=Rc[s+1];
            if (!nl || !nr) continue;
            float cost = C_TRAV + C_ISECT * inv * (L[s]*nl + R[s+1]*nr);
            if (cost < best_cost) { best_cost=cost; best_axis=axis; best_bin=s; best_clo=c_lo; best_scale=scale; }
        }
    }
    if (best_axis < 0) return false;

    // re-partition on winning (axis, bin-boundary)
    float c_lo = 1e30f, c_hi = -1e30f;
    for (auto& p : ps) { float c=p.box.centroid(best_axis); c_lo=std::min(c_lo,c); c_hi=std::max(c_hi,c); }
    float scale = n_bins / (c_hi - c_lo);
    for (auto& p : ps) {
        int b = std::min((int)((p.box.centroid(best_axis) - c_lo) * scale), n_bins - 1);
        (b <= best_bin ? lp : rp).push_back(p);
    }
    return !lp.empty() && !rp.empty();
}

// ─── Binned+A hybrid (threshold T=HYBRID_THRESH) ─────────────────────────────
// N > T : Binned SAH (B=16) 2-way split, recurse with this function
// N ≤ T : Strategy A (4-way cross-axis independent best split)
int build_hybrid_binned_A(std::vector<Node>& ns, const Prims& ps, int d) {
    if ((int)ps.size() <= HYBRID_THRESH) return build_A(ns, ps, d);
    int idx = new_node(ns, ps, d);
    float pa = ns[idx].box.half_perim();
    Prims lp, rp;
    if (!binned_split_2way(ps, pa, 16, lp, rp)) { make_leaf(ns, idx, ps); return idx; }
    int l = build_hybrid_binned_A(ns, lp, d + 1);
    int r = build_hybrid_binned_A(ns, rp, d + 1);
    ns[idx].children = {l, r};
    return idx;
}

// ─── Binned+SAH hybrid (threshold T=HYBRID_THRESH) ───────────────────────────
// N > T : Binned SAH (B=16) 2-way split, recurse with this function
// N ≤ T : 2-way SAH Sweep (precise, exact sweep)
int build_hybrid_binned_sweep(std::vector<Node>& ns, const Prims& ps, int d) {
    if ((int)ps.size() <= HYBRID_THRESH) return build_2way(ns, ps, d);
    int idx = new_node(ns, ps, d);
    float pa = ns[idx].box.half_perim();
    Prims lp, rp;
    if (!binned_split_2way(ps, pa, 16, lp, rp)) { make_leaf(ns, idx, ps); return idx; }
    int l = build_hybrid_binned_sweep(ns, lp, d + 1);
    int r = build_hybrid_binned_sweep(ns, rp, d + 1);
    ns[idx].children = {l, r};
    return idx;
}

// ─── Binned4+A hybrid (threshold T=HYBRID_THRESH) ────────────────────────────
// N > T : two rounds of Binned SAH (B=16) splits inline → up to 4 child groups
//         (no intermediate node recorded — the current node directly owns all groups)
// N ≤ T : Strategy A (4-way cross-axis independent best split)
int build_hybrid_binned4_A(std::vector<Node>& ns, const Prims& ps, int d) {
    if ((int)ps.size() <= HYBRID_THRESH) return build_A(ns, ps, d);
    int idx = new_node(ns, ps, d);
    float pa = ns[idx].box.half_perim();
    // first binned split → lp / rp
    Prims lp, rp;
    if (!binned_split_2way(ps, pa, 16, lp, rp)) { make_leaf(ns, idx, ps); return idx; }
    // second binned split on each half → up to 4 groups total
    std::vector<Prims> groups;
    for (Prims& sub : std::vector<Prims>{lp, rp}) {
        Prims sl, sr;
        float spa = union_box(sub).half_perim();
        if ((int)sub.size() > HYBRID_THRESH && binned_split_2way(sub, spa, 16, sl, sr)) {
            groups.push_back(std::move(sl));
            groups.push_back(std::move(sr));
        } else {
            groups.push_back(std::move(sub));
        }
    }
    for (auto& g : groups) ns[idx].children.push_back(build_hybrid_binned4_A(ns, g, d + 1));
    return idx;
}

// ════════════════════════════════════════════════════════════════════════════
// Collapse k=2  (Fast Quad-BVH baseline)
//   Step 1: build a full 2-way SAH BVH into a temporary array
//   Step 2: walk the temp tree; for every inner node replace its children
//           with its grandchildren, discarding the intermediate layer.
//           If a child is already a leaf it is kept as-is (cannot expand further).
//   Result: ~4-way tree, intermediate nodes gone → lower SAH (no C_trav for them)
// ════════════════════════════════════════════════════════════════════════════
int collapse_k2(const std::vector<Node>& old_ns,
                std::vector<Node>&       new_ns,
                int oi, int d)
{
    int ni = (int)new_ns.size();
    new_ns.push_back(Node());
    new_ns[ni].box   = old_ns[oi].box;
    new_ns[ni].depth = d;

    if (old_ns[oi].is_leaf()) {
        new_ns[ni].prims = old_ns[oi].prims;
        return ni;
    }

    // Collect grandchildren:
    //   child is leaf  → keep it directly (cannot skip into non-existent children)
    //   child is inner → expand one level, exposing its children as grandchildren
    std::vector<int> gc;
    for (int ci : old_ns[oi].children) {
        if (old_ns[ci].is_leaf()) {
            gc.push_back(ci);
        } else {
            for (int gci : old_ns[ci].children) gc.push_back(gci);
        }
    }

    for (int gci : gc) {
        int nc = collapse_k2(old_ns, new_ns, gci, d + 1);
        new_ns[ni].children.push_back(nc);
    }
    return ni;
}

int build_collapse_k2(std::vector<Node>& ns, const Prims& ps, int d) {
    // build 2-way SAH into temporary storage
    std::vector<Node> tmp;
    tmp.reserve(ps.size() * 4);
    build_2way(tmp, ps, 0);
    // collapse: skip every other level
    return collapse_k2(tmp, ns, 0, d);
}

// ─── Binned Collapse k=2 ──────────────────────────────────────────────────────
// Build a 2-way Binned SAH (B=16) tree, then apply level-skipping collapse:
// every inner node adopts its grandchildren, removing the intermediate layer.
// Result: ~4-way tree built from binned splits, same post-process as Collapse k=2.
int build_binned_collapse(std::vector<Node>& ns, const Prims& ps, int d) {
    std::vector<Node> tmp;
    tmp.reserve(ps.size() * 4);
    build_binned(tmp, ps, 0, 16);
    return collapse_k2(tmp, ns, 0, d);
}

// ─── Binned+SAH+Collapse (threshold T=HYBRID_THRESH) ─────────────────────────
// Build the Binned+SAH hybrid tree (2-way throughout), then apply level-skipping
// collapse to the entire tree → ~4-way everywhere, top and bottom included.
int build_hybrid_binned_sweep_collapsed(std::vector<Node>& ns, const Prims& ps, int d) {
    std::vector<Node> tmp;
    tmp.reserve(ps.size() * 4);
    build_hybrid_binned_sweep(tmp, ps, 0);
    return collapse_k2(tmp, ns, 0, d);
}

// ─── Binned4+Coll hybrid (threshold T=HYBRID_THRESH) ─────────────────────────
// Control variable for Binned4+A: identical top (two inline Binned splits → ~4-way)
// but bottom switches to build_collapse_k2 (Sweep SAH + level-skip collapse → ~4-way)
// instead of Strategy A.  Isolates the effect of the bottom-level split method.
int build_hybrid_binned4_coll(std::vector<Node>& ns, const Prims& ps, int d) {
    if ((int)ps.size() <= HYBRID_THRESH) return build_collapse_k2(ns, ps, d);
    int idx = new_node(ns, ps, d);
    float pa = ns[idx].box.half_perim();
    Prims lp, rp;
    if (!binned_split_2way(ps, pa, 16, lp, rp)) { make_leaf(ns, idx, ps); return idx; }
    std::vector<Prims> groups;
    for (Prims& sub : std::vector<Prims>{lp, rp}) {
        Prims sl, sr;
        float spa = union_box(sub).half_perim();
        if ((int)sub.size() > HYBRID_THRESH && binned_split_2way(sub, spa, 16, sl, sr)) {
            groups.push_back(std::move(sl));
            groups.push_back(std::move(sr));
        } else {
            groups.push_back(std::move(sub));
        }
    }
    for (auto& g : groups) ns[idx].children.push_back(build_hybrid_binned4_coll(ns, g, d+1));
    return idx;
}

// ─── JSON helpers ─────────────────────────────────────────────────────────────
static std::string f2s(float v){char buf[32];snprintf(buf,32,"%.6g",v);return buf;}

std::string nodes_to_json(const std::vector<Node>& nodes){
    std::ostringstream o; o<<"[";
    for(int i=0;i<(int)nodes.size();i++){
        const auto& n=nodes[i];
        o<<"{\"id\":"<<i<<",\"depth\":"<<n.depth
         <<",\"box\":["<<f2s(n.box.lo[0])<<","<<f2s(n.box.lo[1])
                  <<","<<f2s(n.box.hi[0])<<","<<f2s(n.box.hi[1])<<"]"
         <<",\"children\":[";
        for(int j=0;j<(int)n.children.size();j++){if(j)o<<",";o<<n.children[j];}
        o<<"],\"prims\":[";
        for(int j=0;j<(int)n.prims.size();j++){if(j)o<<",";o<<n.prims[j];}
        o<<"]}"; if(i+1<(int)nodes.size())o<<",";
    }
    o<<"]"; return o.str();
}

// ─── Main ─────────────────────────────────────────────────────────────────────
int main(){
    struct TestCase { int N; float R; float range; bool run_C; bool is_tri; };
    std::vector<TestCase> tests = {
        { 128,    1.0f,  12.0f, true,  false },  // circle
        { 1024,   1.0f,  36.0f, true,  false },  // circle, C runs (N <= 1024)
        { 102400, 1.0f, 340.0f, false, false },  // circle, C skipped (N > 1024)
        { 128,    1.0f,  12.0f, true,  true  },  // triangle
        { 1024,   1.0f,  36.0f, true,  true  },  // triangle
        { 102400, 1.0f, 340.0f, false, true  },  // triangle
    };

    using BFn = std::function<int(std::vector<Node>&,const Prims&,int)>;
    struct Strategy{ std::string name; BFn fn; bool heavy; };
    std::vector<Strategy> strats={
        {"Median Split",      [](auto& ns,auto& ps,int d){return build_median      (ns,ps,d);},     false},
        {"2-way SAH Sweep",   [](auto& ns,auto& ps,int d){return build_2way        (ns,ps,d);},     false},
        {"Binned SAH (B=16)",      [](auto& ns,auto& ps,int d){return build_binned         (ns,ps,d,16);}, false},
        {"Binned Collapse k=2",    [](auto& ns,auto& ps,int d){return build_binned_collapse(ns,ps,d);},    false},
        {"Collapse k=2",           [](auto& ns,auto& ps,int d){return build_collapse_k2    (ns,ps,d);},    false},
        {"A Independent",   [](auto& ns,auto& ps,int d){return build_A           (ns,ps,d);}, false},
        {"B Hierarchical",  [](auto& ns,auto& ps,int d){return build_B           (ns,ps,d);}, false},
        {"C Exhaustive 2D",   [](auto& ns,auto& ps,int d){return build_C                  (ns,ps,d);},     true },
        {"D All-Axis",        [](auto& ns,auto& ps,int d){return build_D                  (ns,ps,d);},     true },
        {"Binned+A (T=32)",        [](auto& ns,auto& ps,int d){return build_hybrid_binned_A            (ns,ps,d);}, false},
        {"Binned+SAH (T=32)",      [](auto& ns,auto& ps,int d){return build_hybrid_binned_sweep         (ns,ps,d);}, false},
        {"Binned4+A (T=32)",       [](auto& ns,auto& ps,int d){return build_hybrid_binned4_A            (ns,ps,d);}, false},
        {"Binned4+Coll (T=32)",    [](auto& ns,auto& ps,int d){return build_hybrid_binned4_coll          (ns,ps,d);}, false},
        {"Binned+SAH+Coll (T=32)", [](auto& ns,auto& ps,int d){return build_hybrid_binned_sweep_collapsed(ns,ps,d);},false},
    };

    std::mt19937 rng_seed(0);  // master seed

    std::ostringstream json;
    json << "{\"test_cases\":[";

    for(int ti=0;ti<(int)tests.size();ti++){
        auto& tc=tests[ti];
        fprintf(stderr,"\n=== N=%d  R=%.1f  range=[-%.0f,%.0f]  prim=%s ===\n",
                tc.N,tc.R,tc.range,tc.range, tc.is_tri?"triangle":"circle");

        std::mt19937 rng(42 + ti * 1000);
        std::uniform_real_distribution<float> dist(-tc.range, tc.range);
        std::vector<float> cx(tc.N), cy(tc.N);
        // triangle vertex storage: [N][3 vertices][2 coords]
        using TriV = std::array<std::array<float,2>,3>;
        std::vector<TriV> tverts(tc.N);

        Prims base;
        if(tc.is_tri){
            std::uniform_real_distribution<float> vd(-tc.R, tc.R);
            for(int i=0;i<tc.N;i++){
                float ox=dist(rng), oy=dist(rng);
                AABB b;
                for(int v=0;v<3;v++){
                    tverts[i][v][0]=ox+vd(rng);
                    tverts[i][v][1]=oy+vd(rng);
                    b.lo[0]=std::min(b.lo[0],tverts[i][v][0]);
                    b.lo[1]=std::min(b.lo[1],tverts[i][v][1]);
                    b.hi[0]=std::max(b.hi[0],tverts[i][v][0]);
                    b.hi[1]=std::max(b.hi[1],tverts[i][v][1]);
                }
                base.push_back({b, i});
            }
        } else {
            for(int i=0;i<tc.N;i++){cx[i]=dist(rng);cy[i]=dist(rng);}
            for(int i=0;i<tc.N;i++) base.push_back({AABB::circle(cx[i],cy[i],tc.R),i});
        }

        if(ti) json<<",";
        json<<"{\"N\":"<<tc.N<<",\"R\":"<<tc.R<<",\"range\":"<<tc.range;
        if(tc.is_tri){
            json<<",\"type\":\"triangle\",\"triangles\":[";
            for(int i=0;i<tc.N;i++){
                json<<"[["<<f2s(tverts[i][0][0])<<","<<f2s(tverts[i][0][1])<<"],["
                          <<f2s(tverts[i][1][0])<<","<<f2s(tverts[i][1][1])<<"],["
                          <<f2s(tverts[i][2][0])<<","<<f2s(tverts[i][2][1])<<"]]";
                if(i+1<tc.N)json<<",";
            }
            json<<"]";
        } else {
            json<<",\"type\":\"circle\",\"circles\":[";
            for(int i=0;i<tc.N;i++){
                json<<"["<<f2s(cx[i])<<","<<f2s(cy[i])<<"]";
                if(i+1<tc.N)json<<",";
            }
            json<<"]";
        }
        json<<",\"strategies\":[";

        bool first_strat=true;
        for(auto& st:strats){
            if(st.heavy && !tc.run_C){
                fprintf(stderr,"  [%-22s]  SKIPPED (N too large)\n",st.name.c_str());
                continue;
            }
            const int RUNS = 10;
            double dt = 0.0;
            std::vector<Node> nodes;
            for(int run=0;run<RUNS;run++){
                nodes.clear(); nodes.reserve(tc.N*4);
                auto t0=std::chrono::high_resolution_clock::now();
                st.fn(nodes,base,0);
                auto t1=std::chrono::high_resolution_clock::now();
                dt+=std::chrono::duration<double,std::milli>(t1-t0).count();
            }
            dt/=RUNS;
            double cost=sah_cost(nodes);
            auto   s   =tree_stats(nodes);
            assert(s.prim_count==tc.N);

            fprintf(stderr,"  [%-22s]  time=%9.4f ms (avg %d runs)  SAH=%9.4f  nodes=%6d  leaves=%6d  maxdepth=%d\n",
                    st.name.c_str(),dt,RUNS,cost,s.nodes,s.leaves,s.max_depth);

            if(!first_strat) json<<",";
            first_strat=false;
            json<<"{\"name\":\""<<st.name<<"\""
                <<",\"build_time_ms\":"<<dt
                <<",\"sah_cost\":"<<cost
                <<",\"node_count\":"<<s.nodes
                <<",\"leaf_count\":"<<s.leaves
                <<",\"max_depth\":"<<s.max_depth
                <<",\"nodes\":"<<nodes_to_json(nodes)
                <<"}";
        }
        json<<"]}";
    }
    json<<"]}";
    std::cout<<json.str()<<std::endl;
    fprintf(stderr,"\nDone.\n");
    return 0;
}
