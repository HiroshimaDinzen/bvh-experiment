// bvh3d_compute.cpp
// Multi-test BVH comparison: N=128/1024/102400, 3-D version.
// SAH proxy: half surface area  (dx*dy + dy*dz + dz*dx)
// New strategies vs 2D: A4 (4-way cross-axis), A8 (8-way), B4, B8.
// Outputs JSON to stdout, progress to stderr.
// Compile (MSVC): cl /O2 /std:c++17 /EHsc bvh3d_compute.cpp /Fe:bvh3d_compute.exe

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

static const float C_TRAV        = 1.0f;
static const float C_ISECT       = 1.0f;
static const int   MAX_LEAF      = 4;
static const int   HYBRID_THRESH = 32;

// ─── AABB 3-D ─────────────────────────────────────────────────────────────────
struct AABB {
    float lo[3] = {  1e30f,  1e30f,  1e30f };
    float hi[3] = { -1e30f, -1e30f, -1e30f };
    void expand(const AABB& b) {
        for(int i=0;i<3;i++){lo[i]=std::min(lo[i],b.lo[i]);hi[i]=std::max(hi[i],b.hi[i]);}
    }
    // Half surface area: SAH proxy in 3D (equivalent to half-perim in 2D)
    float hsa() const {
        float dx=hi[0]-lo[0], dy=hi[1]-lo[1], dz=hi[2]-lo[2];
        return dx*dy + dy*dz + dz*dx;
    }
    float centroid(int ax) const { return (lo[ax]+hi[ax])*0.5f; }
    static AABB sphere(float cx,float cy,float cz,float r){
        AABB b;
        b.lo[0]=cx-r; b.lo[1]=cy-r; b.lo[2]=cz-r;
        b.hi[0]=cx+r; b.hi[1]=cy+r; b.hi[2]=cz+r;
        return b;
    }
};

struct Prim { AABB box; int idx; };
using  Prims = std::vector<Prim>;

AABB union_box(const Prims& ps){ AABB b; for(auto& p:ps) b.expand(p.box); return b; }

// ─── BVH Node ─────────────────────────────────────────────────────────────────
struct Node {
    AABB             box;
    std::vector<int> children;
    std::vector<int> prims;
    int              depth=0;
    bool is_leaf() const { return !prims.empty(); }
};

// ─── SAH cost (whole tree, uses hsa) ─────────────────────────────────────────
double sah_cost(const std::vector<Node>& nodes){
    float ra=nodes[0].box.hsa();
    if(ra<1e-9f) return 0.0;
    double c=0.0;
    for(auto& n:nodes){
        float w=n.box.hsa()/ra;
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

// ─── 1-D SAH sweep (axis = 0/1/2) ────────────────────────────────────────────
struct Split { int s; float cost; Prims sorted; };

Split best_1d(const Prims& ps, int axis){
    Prims sp=ps;
    std::sort(sp.begin(),sp.end(),[axis](const Prim& a,const Prim& b){
        return a.box.centroid(axis)<b.box.centroid(axis);
    });
    int n=(int)sp.size();
    float pa=union_box(sp).hsa();
    float inv=pa>1e-9f?1.f/pa:0.f;
    std::vector<float> L(n),R(n);
    { AABB b; for(int i=0;i<n;i++)    {b.expand(sp[i].box);L[i]=b.hsa();} }
    { AABB b; for(int i=n-1;i>=0;i--) {b.expand(sp[i].box);R[i]=b.hsa();} }
    float best=C_ISECT*n; int bs=-1;
    for(int s=0;s<n-1;s++){
        float c=C_TRAV+C_ISECT*inv*(L[s]*(s+1)+R[s+1]*(n-s-1));
        if(c<best){best=c;bs=s;}
    }
    return {bs,best,std::move(sp)};
}

// best over all 3 axes
Split best_any(const Prims& ps){
    Split best{-1,C_ISECT*(float)ps.size(),ps};
    for(int ax=0;ax<3;ax++){
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
    float sp[3]={};
    for(int ax=0;ax<3;ax++){
        float lo=1e30f,hi=-1e30f;
        for(auto& p:ps){float c=p.box.centroid(ax);lo=std::min(lo,c);hi=std::max(hi,c);}
        sp[ax]=hi-lo;
    }
    int axis=0; for(int ax=1;ax<3;ax++) if(sp[ax]>sp[axis]) axis=ax;
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

// ─── Binned SAH 3-D (B bins per axis, 3 axes) ────────────────────────────────
int build_binned(std::vector<Node>& ns, const Prims& ps, int d, int n_bins){
    int idx=new_node(ns,ps,d);
    int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    float pa=ns[idx].box.hsa();
    float inv=pa>1e-9f?1.f/pa:0.f;
    float best_cost=C_ISECT*n;
    int   best_axis=-1, best_bin=-1;
    struct Bin{AABB box;int cnt=0;};
    for(int axis=0;axis<3;axis++){
        float c_lo=1e30f,c_hi=-1e30f;
        for(auto& p:ps){float c=p.box.centroid(axis);c_lo=std::min(c_lo,c);c_hi=std::max(c_hi,c);}
        if(c_hi-c_lo<1e-9f) continue;
        float scale=n_bins/(c_hi-c_lo);
        std::vector<Bin> bins(n_bins);
        for(auto& p:ps){int b=std::min((int)((p.box.centroid(axis)-c_lo)*scale),n_bins-1);bins[b].box.expand(p.box);bins[b].cnt++;}
        std::vector<float> L(n_bins); std::vector<int> Lc(n_bins,0);
        std::vector<float> R(n_bins); std::vector<int> Rc(n_bins,0);
        {AABB b;int c=0;for(int i=0;i<n_bins;i++){if(bins[i].cnt>0)b.expand(bins[i].box);c+=bins[i].cnt;L[i]=b.hsa();Lc[i]=c;}}
        {AABB b;int c=0;for(int i=n_bins-1;i>=0;i--){if(bins[i].cnt>0)b.expand(bins[i].box);c+=bins[i].cnt;R[i]=b.hsa();Rc[i]=c;}}
        for(int s=0;s<n_bins-1;s++){
            int nl=Lc[s],nr=Rc[s+1];
            if(!nl||!nr) continue;
            float cost=C_TRAV+C_ISECT*inv*(L[s]*nl+R[s+1]*nr);
            if(cost<best_cost){best_cost=cost;best_axis=axis;best_bin=s;}
        }
    }
    if(best_axis<0){make_leaf(ns,idx,ps);return idx;}
    float c_lo=1e30f,c_hi=-1e30f;
    for(auto& p:ps){float c=p.box.centroid(best_axis);c_lo=std::min(c_lo,c);c_hi=std::max(c_hi,c);}
    float scale=n_bins/(c_hi-c_lo);
    Prims lp,rp;
    for(auto& p:ps){int b=std::min((int)((p.box.centroid(best_axis)-c_lo)*scale),n_bins-1);(b<=best_bin?lp:rp).push_back(p);}
    if(lp.empty()||rp.empty()){make_leaf(ns,idx,ps);return idx;}
    int l=build_binned(ns,lp,d+1,n_bins);
    int r=build_binned(ns,rp,d+1,n_bins);
    ns[idx].children={l,r}; return idx;
}

// ─── A4: 3-D 4-way cross-axis ────────────────────────────────────────────────
// Primary: globally best of 3 axes → L, R
// Secondary: per-half, best of remaining 2 axes → LL, LR, RL, RR (up to 4)
int build_A4(std::vector<Node>& ns, const Prims& ps, int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    // find best axis
    Split sr[3]; for(int ax=0;ax<3;ax++) sr[ax]=best_1d(ps,ax);
    int primary=-1; float best_cost=C_ISECT*n;
    for(int ax=0;ax<3;ax++) if(sr[ax].s>=0&&sr[ax].cost<best_cost){best_cost=sr[ax].cost;primary=ax;}
    if(primary<0){make_leaf(ns,idx,ps);return idx;}
    Prims L=slice(sr[primary].sorted,0,sr[primary].s+1);
    Prims R=slice(sr[primary].sorted,sr[primary].s+1,n);
    // secondary: each half independently picks best of remaining 2 axes
    std::vector<Prims> groups;
    for(Prims* half : {&L,&R}){
        int sec=-1; float sc=C_ISECT*(float)half->size();
        Split best_sec{-1,sc,*half};
        for(int ax=0;ax<3;ax++){
            if(ax==primary) continue;
            auto r=best_1d(*half,ax);
            if(r.s>=0&&r.cost<best_sec.cost){best_sec=std::move(r);sec=ax;}
        }
        if(sec>=0){
            groups.push_back(slice(best_sec.sorted,0,best_sec.s+1));
            groups.push_back(slice(best_sec.sorted,best_sec.s+1,(int)best_sec.sorted.size()));
        } else {
            groups.push_back(*half);
        }
    }
    groups.erase(std::remove_if(groups.begin(),groups.end(),[](const Prims& q){return q.empty();}),groups.end());
    if(groups.size()<=1){make_leaf(ns,idx,ps);return idx;}
    for(auto& g:groups) ns[idx].children.push_back(build_A4(ns,g,d+1));
    return idx;
}

// ─── A8: 3-D 8-way rotating-axis ─────────────────────────────────────────────
// Primary: best of 3 axes → L, R
// Secondary: one axis for both L and R (best combined cost) → LL, LR, RL, RR
// Tertiary: remaining axis applied to each of the 4 groups → up to 8 children
int build_A8(std::vector<Node>& ns, const Prims& ps, int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    // primary
    Split sr[3]; for(int ax=0;ax<3;ax++) sr[ax]=best_1d(ps,ax);
    int primary=-1; float best_cost=C_ISECT*n;
    for(int ax=0;ax<3;ax++) if(sr[ax].s>=0&&sr[ax].cost<best_cost){best_cost=sr[ax].cost;primary=ax;}
    if(primary<0){make_leaf(ns,idx,ps);return idx;}
    Prims L=slice(sr[primary].sorted,0,sr[primary].s+1);
    Prims R=slice(sr[primary].sorted,sr[primary].s+1,n);
    // secondary: pick one axis for both halves (minimize sum of costs on L+R)
    int secondary=-1; float sec_cost=1e30f;
    Split secL[3],secR[3];
    for(int ax=0;ax<3;ax++){
        if(ax==primary) continue;
        secL[ax]=best_1d(L,ax); secR[ax]=best_1d(R,ax);
        float cL=(secL[ax].s>=0)?secL[ax].cost:(C_ISECT*(float)L.size());
        float cR=(secR[ax].s>=0)?secR[ax].cost:(C_ISECT*(float)R.size());
        if(cL+cR<sec_cost){sec_cost=cL+cR;secondary=ax;}
    }
    // tertiary: remaining axis
    int tertiary=-1;
    for(int ax=0;ax<3;ax++) if(ax!=primary&&ax!=secondary){tertiary=ax;break;}
    // apply secondary → up to 4 mid groups
    auto apply=[&](const Split& r, const Prims& half)->std::vector<Prims>{
        if(r.s<0) return {half};
        return {slice(r.sorted,0,r.s+1),slice(r.sorted,r.s+1,(int)r.sorted.size())};
    };
    std::vector<Prims> mids;
    if(secondary<0){mids.push_back(L);mids.push_back(R);}
    else{for(auto& g:apply(secL[secondary],L)) mids.push_back(g);
         for(auto& g:apply(secR[secondary],R)) mids.push_back(g);}
    // apply tertiary → up to 8 groups
    std::vector<Prims> groups;
    if(tertiary<0){groups=mids;}
    else{
        for(auto& m:mids){
            auto rt=best_1d(m,tertiary);
            if(rt.s>=0){groups.push_back(slice(rt.sorted,0,rt.s+1));
                        groups.push_back(slice(rt.sorted,rt.s+1,(int)rt.sorted.size()));}
            else groups.push_back(m);
        }
    }
    groups.erase(std::remove_if(groups.begin(),groups.end(),[](const Prims& q){return q.empty();}),groups.end());
    if(groups.size()<=1){make_leaf(ns,idx,ps);return idx;}
    for(auto& g:groups) ns[idx].children.push_back(build_A8(ns,g,d+1));
    return idx;
}

// ─── B4: 3-D hierarchical 4-way (2 greedy levels, any axis) ──────────────────
int build_B4(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    auto r1=best_any(ps);
    if(r1.s<0){make_leaf(ns,idx,ps);return idx;}
    Prims L=slice(r1.sorted,0,r1.s+1);
    Prims R=slice(r1.sorted,r1.s+1,n);
    auto r2=best_any(L),r3=best_any(R);
    std::vector<Prims> quads;
    if(r2.s>=0){quads.push_back(slice(r2.sorted,0,r2.s+1));
                quads.push_back(slice(r2.sorted,r2.s+1,(int)r2.sorted.size()));}
    else quads.push_back(L);
    if(r3.s>=0){quads.push_back(slice(r3.sorted,0,r3.s+1));
                quads.push_back(slice(r3.sorted,r3.s+1,(int)r3.sorted.size()));}
    else quads.push_back(R);
    quads.erase(std::remove_if(quads.begin(),quads.end(),[](const Prims& q){return q.empty();}),quads.end());
    if(quads.size()<=1){make_leaf(ns,idx,ps);return idx;}
    for(auto& q:quads) ns[idx].children.push_back(build_B4(ns,q,d+1));
    return idx;
}

// ─── B8: 3-D hierarchical 8-way (3 greedy levels, any axis) ──────────────────
int build_B8(std::vector<Node>& ns,const Prims& ps,int d){
    int idx=new_node(ns,ps,d); int n=(int)ps.size();
    if(n<=MAX_LEAF){make_leaf(ns,idx,ps);return idx;}
    auto r1=best_any(ps);
    if(r1.s<0){make_leaf(ns,idx,ps);return idx;}
    Prims L=slice(r1.sorted,0,r1.s+1);
    Prims R=slice(r1.sorted,r1.s+1,n);
    auto r2=best_any(L),r3=best_any(R);
    std::vector<Prims> mids;
    if(r2.s>=0){mids.push_back(slice(r2.sorted,0,r2.s+1));
                mids.push_back(slice(r2.sorted,r2.s+1,(int)r2.sorted.size()));}
    else mids.push_back(L);
    if(r3.s>=0){mids.push_back(slice(r3.sorted,0,r3.s+1));
                mids.push_back(slice(r3.sorted,r3.s+1,(int)r3.sorted.size()));}
    else mids.push_back(R);
    std::vector<Prims> octs;
    for(auto& m:mids){
        auto rm=best_any(m);
        if(rm.s>=0){octs.push_back(slice(rm.sorted,0,rm.s+1));
                    octs.push_back(slice(rm.sorted,rm.s+1,(int)rm.sorted.size()));}
        else octs.push_back(m);
    }
    octs.erase(std::remove_if(octs.begin(),octs.end(),[](const Prims& q){return q.empty();}),octs.end());
    if(octs.size()<=1){make_leaf(ns,idx,ps);return idx;}
    for(auto& g:octs) ns[idx].children.push_back(build_B8(ns,g,d+1));
    return idx;
}

// ─── Hybrid helpers (3-axis version) ─────────────────────────────────────────
static bool binned_split_2way(const Prims& ps, float pa, int n_bins, Prims& lp, Prims& rp){
    int n=(int)ps.size();
    float inv=pa>1e-9f?1.f/pa:0.f;
    float best_cost=C_ISECT*n;
    int   best_axis=-1,best_bin=-1;
    float best_clo=0.f,best_scale=0.f;
    struct Bin{AABB box;int cnt=0;};
    for(int axis=0;axis<3;axis++){
        float c_lo=1e30f,c_hi=-1e30f;
        for(auto& p:ps){float c=p.box.centroid(axis);c_lo=std::min(c_lo,c);c_hi=std::max(c_hi,c);}
        if(c_hi-c_lo<1e-9f) continue;
        float scale=n_bins/(c_hi-c_lo);
        std::vector<Bin> bins(n_bins);
        for(auto& p:ps){int b=std::min((int)((p.box.centroid(axis)-c_lo)*scale),n_bins-1);bins[b].box.expand(p.box);bins[b].cnt++;}
        std::vector<float> L(n_bins);std::vector<int> Lc(n_bins,0);
        std::vector<float> R(n_bins);std::vector<int> Rc(n_bins,0);
        {AABB b;int c=0;for(int i=0;i<n_bins;i++){if(bins[i].cnt>0)b.expand(bins[i].box);c+=bins[i].cnt;L[i]=b.hsa();Lc[i]=c;}}
        {AABB b;int c=0;for(int i=n_bins-1;i>=0;i--){if(bins[i].cnt>0)b.expand(bins[i].box);c+=bins[i].cnt;R[i]=b.hsa();Rc[i]=c;}}
        for(int s=0;s<n_bins-1;s++){
            int nl=Lc[s],nr=Rc[s+1];
            if(!nl||!nr) continue;
            float cost=C_TRAV+C_ISECT*inv*(L[s]*nl+R[s+1]*nr);
            if(cost<best_cost){best_cost=cost;best_axis=axis;best_bin=s;best_clo=c_lo;best_scale=scale;}
        }
    }
    if(best_axis<0) return false;
    float c_lo=1e30f,c_hi=-1e30f;
    for(auto& p:ps){float c=p.box.centroid(best_axis);c_lo=std::min(c_lo,c);c_hi=std::max(c_hi,c);}
    float scale=n_bins/(c_hi-c_lo);
    for(auto& p:ps){int b=std::min((int)((p.box.centroid(best_axis)-c_lo)*scale),n_bins-1);(b<=best_bin?lp:rp).push_back(p);}
    return !lp.empty()&&!rp.empty();
}

// Binned+A4 (T=32): top Binned 2-way, bottom A4
int build_hybrid_binned_A4(std::vector<Node>& ns,const Prims& ps,int d){
    if((int)ps.size()<=HYBRID_THRESH) return build_A4(ns,ps,d);
    int idx=new_node(ns,ps,d);
    float pa=ns[idx].box.hsa();
    Prims lp,rp;
    if(!binned_split_2way(ps,pa,16,lp,rp)){make_leaf(ns,idx,ps);return idx;}
    int l=build_hybrid_binned_A4(ns,lp,d+1);
    int r=build_hybrid_binned_A4(ns,rp,d+1);
    ns[idx].children={l,r}; return idx;
}

// Binned4+A4 (T=32): top 2-round inline Binned (~4-way), bottom A4
int build_hybrid_binned4_A4(std::vector<Node>& ns,const Prims& ps,int d){
    if((int)ps.size()<=HYBRID_THRESH) return build_A4(ns,ps,d);
    int idx=new_node(ns,ps,d);
    float pa=ns[idx].box.hsa();
    Prims lp,rp;
    if(!binned_split_2way(ps,pa,16,lp,rp)){make_leaf(ns,idx,ps);return idx;}
    std::vector<Prims> groups;
    for(Prims& sub:std::vector<Prims>{lp,rp}){
        Prims sl,sr;
        float spa=union_box(sub).hsa();
        if((int)sub.size()>HYBRID_THRESH&&binned_split_2way(sub,spa,16,sl,sr)){
            groups.push_back(std::move(sl));
            groups.push_back(std::move(sr));
        } else {
            groups.push_back(std::move(sub));
        }
    }
    for(auto& g:groups) ns[idx].children.push_back(build_hybrid_binned4_A4(ns,g,d+1));
    return idx;
}

// Binned4+A8 (T=32): top 2-round inline Binned (~4-way), bottom A8
int build_hybrid_binned4_A8(std::vector<Node>& ns,const Prims& ps,int d){
    if((int)ps.size()<=HYBRID_THRESH) return build_A8(ns,ps,d);
    int idx=new_node(ns,ps,d);
    float pa=ns[idx].box.hsa();
    Prims lp,rp;
    if(!binned_split_2way(ps,pa,16,lp,rp)){make_leaf(ns,idx,ps);return idx;}
    std::vector<Prims> groups;
    for(Prims& sub:std::vector<Prims>{lp,rp}){
        Prims sl,sr;
        float spa=union_box(sub).hsa();
        if((int)sub.size()>HYBRID_THRESH&&binned_split_2way(sub,spa,16,sl,sr)){
            groups.push_back(std::move(sl));
            groups.push_back(std::move(sr));
        } else {
            groups.push_back(std::move(sub));
        }
    }
    for(auto& g:groups) ns[idx].children.push_back(build_hybrid_binned4_A8(ns,g,d+1));
    return idx;
}

// ════════════════════════════════════════════════════════════════════════════
// Collapse k=2
// ════════════════════════════════════════════════════════════════════════════
int collapse_k2(const std::vector<Node>& old_ns,std::vector<Node>& new_ns,int oi,int d){
    int ni=(int)new_ns.size();
    new_ns.push_back(Node());
    new_ns[ni].box=old_ns[oi].box;
    new_ns[ni].depth=d;
    if(old_ns[oi].is_leaf()){new_ns[ni].prims=old_ns[oi].prims;return ni;}
    std::vector<int> gc;
    for(int ci:old_ns[oi].children){
        if(old_ns[ci].is_leaf()) gc.push_back(ci);
        else for(int gci:old_ns[ci].children) gc.push_back(gci);
    }
    for(int gci:gc){int nc=collapse_k2(old_ns,new_ns,gci,d+1);new_ns[ni].children.push_back(nc);}
    return ni;
}

int build_collapse_k2(std::vector<Node>& ns,const Prims& ps,int d){
    std::vector<Node> tmp; tmp.reserve(ps.size()*4);
    build_2way(tmp,ps,0);
    return collapse_k2(tmp,ns,0,d);
}

int build_binned_collapse(std::vector<Node>& ns,const Prims& ps,int d){
    std::vector<Node> tmp; tmp.reserve(ps.size()*4);
    build_binned(tmp,ps,0,16);
    return collapse_k2(tmp,ns,0,d);
}

// ─── JSON helpers ─────────────────────────────────────────────────────────────
static std::string f2s(float v){char buf[32];snprintf(buf,32,"%.6g",v);return buf;}

std::string nodes_to_json(const std::vector<Node>& nodes){
    std::ostringstream o; o<<"[";
    for(int i=0;i<(int)nodes.size();i++){
        const auto& n=nodes[i];
        o<<"{\"id\":"<<i<<",\"depth\":"<<n.depth
         <<",\"box\":["<<f2s(n.box.lo[0])<<","<<f2s(n.box.lo[1])<<","<<f2s(n.box.lo[2])
                  <<","<<f2s(n.box.hi[0])<<","<<f2s(n.box.hi[1])<<","<<f2s(n.box.hi[2])<<"]"
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
    struct TestCase{int N;float R;float range;bool is_tri;};
    std::vector<TestCase> tests={
        {128,   1.0f, 12.0f, false},
        {1024,  1.0f, 36.0f, false},
        {102400,1.0f,340.0f, false},
        {128,   1.0f, 12.0f, true },
        {1024,  1.0f, 36.0f, true },
        {102400,1.0f,340.0f, true },
    };

    using BFn=std::function<int(std::vector<Node>&,const Prims&,int)>;
    struct Strategy{std::string name;BFn fn;};
    std::vector<Strategy> strats={
        {"Median Split",        [](auto& ns,auto& ps,int d){return build_median          (ns,ps,d);}},
        {"2-way SAH Sweep",     [](auto& ns,auto& ps,int d){return build_2way            (ns,ps,d);}},
        {"Binned SAH (B=16)",   [](auto& ns,auto& ps,int d){return build_binned          (ns,ps,d,16);}},
        {"Collapse k=2",        [](auto& ns,auto& ps,int d){return build_collapse_k2     (ns,ps,d);}},
        {"Binned Collapse k=2", [](auto& ns,auto& ps,int d){return build_binned_collapse (ns,ps,d);}},
        {"A4 Independent",      [](auto& ns,auto& ps,int d){return build_A4              (ns,ps,d);}},
        {"A8 Independent",      [](auto& ns,auto& ps,int d){return build_A8              (ns,ps,d);}},
        {"B4 Hierarchical",     [](auto& ns,auto& ps,int d){return build_B4              (ns,ps,d);}},
        {"B8 Hierarchical",     [](auto& ns,auto& ps,int d){return build_B8              (ns,ps,d);}},
        {"Binned+A4 (T=32)",    [](auto& ns,auto& ps,int d){return build_hybrid_binned_A4 (ns,ps,d);}},
        {"Binned4+A4 (T=32)",   [](auto& ns,auto& ps,int d){return build_hybrid_binned4_A4(ns,ps,d);}},
        {"Binned4+A8 (T=32)",   [](auto& ns,auto& ps,int d){return build_hybrid_binned4_A8(ns,ps,d);}},
    };

    std::ostringstream json;
    json<<"{\"test_cases\":[";

    for(int ti=0;ti<(int)tests.size();ti++){
        auto& tc=tests[ti];
        fprintf(stderr,"\n=== N=%d  R=%.1f  range=[-%.0f,%.0f]  prim=%s ===\n",
                tc.N,tc.R,tc.range,tc.range,tc.is_tri?"triangle3d":"sphere");

        std::mt19937 rng(42+ti*1000);
        std::uniform_real_distribution<float> dist(-tc.range,tc.range);

        using TriV=std::array<std::array<float,3>,3>;
        std::vector<float> cx(tc.N),cy(tc.N),cz(tc.N);
        std::vector<TriV> tverts(tc.N);

        Prims base;
        if(tc.is_tri){
            std::uniform_real_distribution<float> vd(-tc.R,tc.R);
            for(int i=0;i<tc.N;i++){
                float ox=dist(rng),oy=dist(rng),oz=dist(rng);
                AABB b;
                for(int v=0;v<3;v++){
                    tverts[i][v][0]=ox+vd(rng);
                    tverts[i][v][1]=oy+vd(rng);
                    tverts[i][v][2]=oz+vd(rng);
                    for(int k=0;k<3;k++){
                        b.lo[k]=std::min(b.lo[k],tverts[i][v][k]);
                        b.hi[k]=std::max(b.hi[k],tverts[i][v][k]);
                    }
                }
                base.push_back({b,i});
            }
        } else {
            for(int i=0;i<tc.N;i++){cx[i]=dist(rng);cy[i]=dist(rng);cz[i]=dist(rng);}
            for(int i=0;i<tc.N;i++) base.push_back({AABB::sphere(cx[i],cy[i],cz[i],tc.R),i});
        }

        if(ti) json<<",";
        json<<"{\"N\":"<<tc.N<<",\"R\":"<<tc.R<<",\"range\":"<<tc.range;
        if(tc.is_tri){
            json<<",\"type\":\"triangle3d\",\"triangles\":[";
            for(int i=0;i<tc.N;i++){
                json<<"[["<<f2s(tverts[i][0][0])<<","<<f2s(tverts[i][0][1])<<","<<f2s(tverts[i][0][2])<<"],["
                          <<f2s(tverts[i][1][0])<<","<<f2s(tverts[i][1][1])<<","<<f2s(tverts[i][1][2])<<"],["
                          <<f2s(tverts[i][2][0])<<","<<f2s(tverts[i][2][1])<<","<<f2s(tverts[i][2][2])<<"]]";
                if(i+1<tc.N) json<<",";
            }
            json<<"]";
        } else {
            json<<",\"type\":\"sphere\",\"spheres\":[";
            for(int i=0;i<tc.N;i++){
                json<<"["<<f2s(cx[i])<<","<<f2s(cy[i])<<","<<f2s(cz[i])<<"]";
                if(i+1<tc.N) json<<",";
            }
            json<<"]";
        }
        json<<",\"strategies\":[";

        bool first_strat=true;
        for(auto& st:strats){
            const int RUNS=10;
            double dt=0.0;
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
