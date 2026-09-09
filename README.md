# BVH Construction Strategy Comparison

Comparing BVH construction strategies in 2D and 3D, including a proposed
method that builds a wide (4-ary / 8-ary) node **directly** during
construction instead of building a binary tree and collapsing it afterwards.

---

## Read this first: what is and is not claimed

### The claim

At **matched arity and a stated SIMD width**, a hybrid of an inline-widened
binned top and a multi-axis grid-partitioned bottom (`Binned8+Ours8`,
`Binned4+Ours4`) builds **faster** than the equivalent collapse-based
baseline, with SAH cost **better in 2D (4.6–10 %)** and **roughly equal in
3D (−0.8 % to +2.2 %)**.

Build-time evidence, stated by strength:

- **3D: solid.** −15.2 % on both meshes. Dragon 630.1 ± 7.5 vs 742.8 ± 7.3 ms
  (10.7 SD apart); Bunny 92.2 ± 1.9 vs 108.7 ± 1.0 ms (7.7 SD).
- **2D: directionally consistent but noisy.** All four configurations favour
  the proposal (−14.8 % to −23.0 %), but two of them separate by only
  ~1.5–1.7 SD. Run-to-run variation between invocations of the same binary
  reached 7 %, which is larger than some of the effects being claimed. The 2D
  build-time advantage should be treated as suggestive, not established.

SAH figures are deterministic and reproduced bit-for-bit across runs; only
the timing carries this uncertainty.

### Scope: the proposal only pays off at the bottom of a hybrid

`Ours4` / `Ours8` are **not** competitive as standalone builders, and this is
not a tuning problem — it is structural.

Per node, `Ours` does the same work as an exact sweep: three sorts and
`3(n−1)` candidate evaluations. Binned SAH does `O(n + B)` with 45
evaluations and **no sort at all**. So `Ours` reduces the *number* of split
levels (Dragon: max depth 25 → 15 for `Ours4`) while raising the cost of
*each* level by more than the depth saving returns:

| Method | Per-node work | Evaluations | Max depth | Build [ms] |
|---|---|---|---|---|
| 2-way Sweep | 3 sorts, `O(n log n)` | `3(n−1)` | 24 | 2115 |
| Binned (B=16) | binning, `O(n+B)`, no sort | 45 | 25 | **682** |
| Ours4 | 3 sorts, `O(n log n)` | `3(n−1)` | **15** | 1370 |
| Ours8 | 3 sorts, `O(n log n)` | `3(n−1)` | **11** | 1029 |

`Ours4` cuts depth by 40 % and is still **2× slower than binned**. Counting
split rounds without pricing each round is exactly the kind of cost
evaluation this repository got wrong before.

The depth argument holds only where the sort is cheap — at the bottom of a
hybrid, where `n ≤ T = 32`. That is where the measured gain comes from:
`Binned8 inline` (binned all the way down) 722.7 ms → `Binned8+Ours8` (Ours
below T) **630.1 ms**, a 13 % saving. Applied to a whole tree the same
argument fails.

**So the claim is scoped to: `Ours` as the bottom stage of a hybrid whose top
is an inline-widened binned build of matching arity.** Nothing broader is
claimed, and the standalone numbers above should not be read as a fallback
position.

### The retraction

An earlier version of this work claimed a **"40–60 % SAH reduction."**
**That claim is withdrawn and was wrong.** It came from comparing an 8-ary
proposal against a 4-ary baseline under a cost model that charges one
traversal cost per inner node regardless of how many children it has. Nearly
all of that apparent gain was the effect of *widening the tree*, which the
existing collapse method (Dammertz et al.) also achieves — and achieves
slightly better. Any document in this repository still asserting a 40–60 %
reduction is superseded by this file.

### What is explicitly not shown

- **No real ray traversal time was measured.** Every SAH number here is a
  proxy metric. The build-time advantage is measured wall-clock time; the
  SAH comparison is a model.
- **The standalone proposals (`Ours4`, `Ours8`) are dominated** by the
  collapse baseline — worse SAH *and* slower to build. Only the hybrids
  survive scrutiny.
- **The finished tree is larger**: 26 % more nodes in 2D, 47 % in 3D. Peak
  memory *during* construction is roughly 40 % lower (no intermediate binary
  tree is materialised), but the resulting tree costs more memory to store.
- **2D and 3D disagree on SAH** (2D consistently better, 3D roughly level).
  This difference is not currently explained.

---

## Assumptions and preconditions

Every number in this repository depends on the following. They are stated
here because the central error in the earlier version of this work was
leaving the first one unstated.

### 1. SIMD width — the load-bearing assumption

The whole-tree SAH cost is

```
SAH = C_trav · Σ_{n ∈ Inner} M(n)/M(root)  +  C_isect · Σ_{n ∈ Leaf} N(n)·M(n)/M(root)
```

where `M` is the half-perimeter (2D) or half surface area (3D) of a node's
bounding box, and `N(n)` is the primitive count of leaf `n`.

This form charges **one** `C_trav` per inner node **regardless of its
arity** `k`. That is only correct when the SIMD width `W` is at least `k`,
so that testing `k` child boxes is a single vector operation. It is the
standard model in the wide-BVH literature (Dammertz et al. 2008; Ernst &
Greiner 2008) — but it is a model *with a precondition*, not a neutral
metric.

**Consequence: SAH values from trees of different arity are not comparable
unless a width `W` is stated.** `sah_cost_simd(nodes, W)` charges
`ceil(k/W)` vector box tests per `k`-ary node and is the honest comparison.
Setting `W ≥ max arity` recovers the flat model above.

The benchmarks print `maxK` (the largest arity actually produced) so that
`W ≥ k` can be **checked**, not assumed. All reported results satisfy it:
4-ary methods reach `maxK = 4` and are compared at `W = 4`; 8-ary methods
reach `maxK = 8` and are compared at `W = 8`.

At the wrong width the ranking inverts. On Dragon, `Binned8+Ours8` scores
24.32 at `W = 8` but 36.07 at `W = 4`, where the 4-ary
`Binned Collapse k=2` (30.48) beats it. **An 8-ary method has no advantage
on 4-wide hardware.**

### 2. Fixed parameters

| Parameter | Value |
|---|---|
| `C_TRAV`, `C_ISECT` | 1.0, 1.0 |
| `MAX_LEAF` | 4 primitives |
| Bin count `B` | 16 |
| Hybrid threshold `T` | 32 (runtime-settable with `--T=<n>`) |
| Area proxy | half-perimeter (2D); half surface area (3D) |

### 3. Timing method

Untimed warm-up builds run first, then `n` timed builds; the mean and
**sample standard deviation** are reported. 3D mesh runs use `n = 7`
(`--runs=<n>`), 2D uses `n = 10`. Differences smaller than a few SD are not
treated as real. Timings are single-threaded and machine-specific; only
*relative* differences measured in the same run are meaningful.

### 4. Arity must be matched before comparing

A hybrid whose top is 4-ary but whose bottom is 8-ary is not a fair
opponent for a uniformly 8-ary tree — SAH weight concentrates near the root,
so a narrow top wastes lanes exactly where cost is highest. Fixing this
(`Binned4+Ours8` → `Binned8+Ours8`, i.e. two inline binned rounds → three)
moved Dragon's SAH from 29.60 to 24.32. Most of the previously reported gap
was this mismatch, not the method.

---

## The proposal

**Multi-axis grid partitioning.** At one node, run a single 1-D SAH sweep
per coordinate axis to obtain one global threshold per axis, then make a
single linear pass that turns the per-axis comparisons into a bit key naming
the child. This yields a 2×2 (`Ours4`) or 2×2×2 (`Ours8`) grid partition.
The sweeps happen once per node and are never repeated for the children.

**No binary tree is built and no collapse pass is run** — the wide node is
committed to directly. This is the structural difference from the wide BVHs
of Dammertz et al., Wald et al., and Ylitie et al., all of which derive
their wide tree *after* a complete binary construction.

**The hybrids** use an inline-widened binned top down to `N ≤ T`, then
switch to `Ours` beneath. The top performs `R` successive binned rounds
without recording the intermediate nodes, producing a `2^R`-ary node:
`R = 2` gives a 4-ary top, `R = 3` an 8-ary top.

### Verified: inline widening ≡ collapsing

`R` inline binned rounds produce a tree **topologically identical** to
building a binary binned tree and collapsing it by `R−1` levels. SAH, node
count, leaf count and max depth match **exactly** — verified for `R = 2` and
`R = 3`, in 2D and 3D, across every test configuration.

This matters for the experiment design: it makes `Binned8+Ours8` versus
`Binned Collapse 8ary` a **clean ablation** with an identical top, differing
only in the bottom. It also isolates where the build-time saving comes from:
the intermediate binary tree is never materialised.

---

## Results at matched arity

### 3D, real meshes, `W = 8`, `T = 32`, n = 7

Dragon (871,306 triangles):

| Method | SAH (W=8) | Build [ms] | Nodes | maxK |
|---|---|---|---|---|
| Binned SAH (B=16), binary baseline | 50.163 | 681.6 ± 6.8 | 575,115 | 2 |
| Binned Collapse 8ary (existing) | **23.811** | 742.8 ± 7.3 | 385,715 | 8 |
| Binned8 inline, pure (same topology) | 23.811 | 722.7 ± 5.7 | 385,715 | 8 |
| **Binned8+Ours8 (proposed)** | 24.323 (+2.2 %) | **630.1 ± 7.5 (−15.2 %)** | 565,372 | 8 |
| Ours8 standalone | 26.633 | 1028.7 ± 7.4 | 553,901 | 8 |

Bunny (144,046 triangles):

| Method | SAH (W=8) | Build [ms] | Nodes | maxK |
|---|---|---|---|---|
| Binned SAH (B=16), binary baseline | 37.759 | 101.6 ± 0.9 | 91,773 | 2 |
| Binned Collapse 8ary (existing) | 18.714 | 108.7 ± 1.0 | 62,865 | 8 |
| **Binned8+Ours8 (proposed)** | **18.560 (−0.8 %)** | **92.2 ± 1.9 (−15.2 %)** | 95,074 | 8 |
| Ours8 standalone | 19.220 | 150.2 ± 3.1 | 91,101 | 8 |

### 2D, synthetic, `W = 4`, `T = 32`, n = 10

`Binned4+Ours4` versus `Binned Collapse k=2`:

| Configuration | ΔSAH | ΔBuild time | Timing separation |
|---|---|---|---|
| N=1024, circles | −4.6 % | −18.2 % | ~1.5 SD (marginal) |
| N=102400, circles | −7.7 % | −17.2 % | ~3.3 SD |
| N=1024, triangles | −5.7 % | −23.0 % | ~13 SD |
| N=102400, triangles | −10.0 % | −14.8 % | ~1.7 SD (marginal) |

Negative is better. In 2D the proposal wins on **both** axes in every
configuration, but see the caveat above: the ΔSAH column is exact and
reproducible, while two of the four timing differences are within ~1.7 SD.

---

## Files

| File | Contents |
|---|---|
| `bvh_compute.cpp` | 2D benchmark, all strategies, emits JSON |
| `bvh3d_compute.cpp` | 3D benchmark, synthetic + OBJ meshes |
| `bvh_compare.py`, `bvh3d_compare.py`, `bvh4_compare.py` | plotting / comparison |
| `bvh_sweep_vis.py` | sweep visualisation |
| `bvh_strategies.md` (`_en`, `_ja`) | full per-strategy notes and raw logs, 3 languages |

Key functions in `bvh3d_compute.cpp`:

- `sah_cost()` — flat cost; valid only for `W ≥ maxK`
- `sah_cost_simd(nodes, W)` — width-aware cost, `ceil(k/W)` per node
- `collapse_n(..., lvl)` — generalised collapse to `2^(lvl+1)`-ary
- `build_binnedN_inline(..., R)` — pure `R`-round inline binned (the control)
- `build_hybrid_binnedN_ours8(..., R)` — the proposed hybrid

## Building and running

```bash
cl /O2 /std:c++17 /EHsc /Fe:bvh3d_compute.exe bvh3d_compute.cpp
bvh3d_compute.exe meshes/dragon/dragon.obj --T=32 --runs=7
```

```bash
cl /O2 /std:c++17 /EHsc /Fe:bvh_compute.exe bvh_compute.cpp
bvh_compute.exe --T=32
```

## Open problems

1. **Measure real ray traversal time.** SAH is a proxy; a 2 % SAH difference
   should not be argued from the model alone.
2. **Explain the 2D/3D discrepancy** in SAH direction.
3. **Make the split decision match the structure built.** Thresholds are
   currently chosen by per-axis *binary* SAH, but a 4-way or 8-way node is
   then built; the cost of that wide split is never evaluated against the
   binary alternative. The criterion and the constructed structure do not
   correspond.
4. **Account for the larger finished tree** (26–47 % more nodes) in any
   memory argument.
