# BVH Construction Strategy Comparison

Comparing BVH construction strategies in 2D and 3D, including a proposed
method that builds a wide (4-ary / 8-ary) node **directly** during
construction instead of building a binary tree and collapsing it afterwards.

---

## Read this first: what is and is not claimed

### The claim

Building a wide BVH **directly**, with no intermediate binary tree, is a
**build-time / quality trade — not a quality win**. At matched arity, a
stated SIMD width, and each builder at *its own* optimal leaf size,
`Binned8+Ours8` gives up **4–14 % SAH** and builds **5–31 % faster** than
the equivalent collapse baseline.

It is also **markedly less sensitive to the leaf-size parameter**: over
`MAX_LEAF` 1→16 its SAH moves +46 % and its build time −35 %, against
+117 % and −54 % for collapse.

**The proposal does not have lower SAH cost.** Swept over `C_trav:C_isect`
of 1:8 … 4:1 and `MAX_LEAF` of 1…16, on two meshes, collapse wins on SAH at
every tuned configuration. See "Parameter sweeps" below — this was tested
hard and the answer did not change.

Both meshes, each method at its own optimum (`MAX_LEAF = 1`):

| `C_trav` | Dragon ΔSAH / Δbuild | Bunny ΔSAH / Δbuild |
|---|---|---|
| 0.125 | +6.2 % / **−24.8 %** | +7.5 % / **−31.0 %** |
| 0.5 | — | +14.1 % / −31.1 % |
| 1 | +8.8 % / −24.2 % | +9.9 % / −29.4 % |
| 2 | — | +8.3 % / −21.3 % |
| 4 | +3.9 % / −7.2 % | +6.2 % / −4.7 % |

Positive ΔSAH means the proposal is worse. The trade is best where traversal
is cheap relative to intersection and worst where it is expensive — at
`C_trav = 4` the build-time advantage nearly vanishes.

**Why the build-time advantage is largest where quality wants small leaves:**
collapse must materialise the entire binary tree before merging levels away,
so the deeper that tree, the more levels `Ours` skips by splitting 8 ways
directly.

Evidence quality: SAH figures are deterministic and reproduce bit-for-bit.
3D build times are solid (Dragon 10.7 SD, Bunny 7.7 SD at `MAX_LEAF = 4`).
2D build times are directionally consistent across all four configurations
but two separate by only ~1.5–1.7 SD, and run-to-run drift between
invocations of the same binary reached 7 % — treat 2D timing as suggestive,
not established.

### Superseded: the earlier "SAH roughly equal in 3D" reading

A previous version of this file reported SAH as roughly level in 3D (−0.8 %
to +2.2 %). That was measured at `MAX_LEAF = 4`, which the sweep later showed
sits almost exactly on the crossover point where the sign flips. It was a
coincidence of an arbitrary parameter, not a property of the method. At the
SAH-optimal leaf size the gap is a consistent 4–14 % in collapse's favour.

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
- **The 2D results have not been re-swept.** 2D still reports SAH better by
  4.6–10 % at `MAX_LEAF = 4`, but that is the exact parameter value the 3D
  sweep showed to be a crossover artifact. Until 2D is swept over `MAX_LEAF`
  and `C_trav` the same way, treat the 2D SAH advantage as unconfirmed.
- **Only two meshes.** Bunny and Dragon. They agree on the sweep results,
  but two meshes is not a survey.

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
| `C_TRAV`, `C_ISECT` | 1.0, 1.0 default; `C_TRAV` settable with `--ctrav=<f>` |
| `MAX_LEAF` | 4 default; settable with `--leaf=<n>` |
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

## Parameter sweeps

Two parameters that had been fixed by convention turned out to control the
verdict, so both were swept. `C_ISECT` is held at 1; `C_trav` is the ratio.

### `MAX_LEAF`, at `C_trav = 1` (SAH at `W = 8`, Ours vs collapse)

| `MAX_LEAF` | Dragon ΔSAH / Δbuild | Bunny ΔSAH / Δbuild |
|---|---|---|
| 1 | +8.8 % / −22.8 % | +9.9 % / −29.0 % |
| 2 | +8.5 % / −24.1 % | +9.6 % / −25.1 % |
| 4 | +2.2 % / −15.3 % | −0.8 % / −14.0 % |
| 8 | −7.1 % / −7.8 % | −6.4 % / −10.3 % |
| 16 | −26.5 % / +10.3 % | −30.0 % / +13.3 % |

The SAH gap is monotone in `MAX_LEAF` and changes sign near 4. **`MAX_LEAF`
must not be raised to make the proposal look better**: doing so degrades both
trees in absolute terms (Dragon collapse 23.811 → 46.562) and the proposal
only "wins" because collapse degrades faster. Collapse inherits the binary
tree's leaves and fills to the cap (11.22 primitives per leaf at 16); Ours
self-regulates through its own SAH termination (4.76).

### `C_trav` × `MAX_LEAF` on Bunny (ΔSAH, positive = proposal worse)

| `C_trav` \ `MAX_LEAF` | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 0.125 | **+7.5** | −3.7 | −16.9 | −17.0 | −41.6 |
| 0.5 | **+14.1** | +5.1 | −7.5 | −11.4 | −36.0 |
| 1 | **+9.9** | +9.6 | −0.8 | −6.4 | −30.0 |
| 2 | **+8.3** | +8.3 | +6.5 | −0.4 | −21.3 |
| 4 | **+6.2** | +6.2 | +6.2 | +5.4 | −11.0 |

Bold is each row's optimal leaf size for both builders, which is 1 at every
ratio tested. Every cell where the proposal wins lies at a leaf size that is
suboptimal for both.

At `C_trav = 4` the SAH termination criterion fires before the leaf cap is
reached, so `MAX_LEAF` 1/2/4 give identical trees and the parameter stops
mattering at all — and collapse still wins by 6.2 %. **The leaf cap cannot
explain the gap.**

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

1. **Measure real ray traversal time.** SAH is a proxy for node visits and
   primitive tests; measuring those directly tests the model instead of
   assuming it. Note the finished tree carries 26–47 % more nodes, which SAH
   discounts (they sit in small boxes) but cache pressure does not — expect
   measured traversal to be worse than SAH predicts.
2. **Explain the 2D/3D discrepancy** in SAH direction. Note the 3D side of
   that comparison was measured at `MAX_LEAF = 4`; re-check it at the swept
   optimum before treating the discrepancy as real.
3. **Make the split decision match the structure built.** Thresholds are
   currently chosen by per-axis *binary* SAH, but a 4-way or 8-way node is
   then built; the cost of that wide split is never evaluated against the
   binary alternative. The criterion and the constructed structure do not
   correspond.
4. **Account for the larger finished tree** (26–47 % more nodes) in any
   memory argument.
