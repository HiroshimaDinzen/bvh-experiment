# BVH Split Strategy Comparison Experiment

> **Read this section first.** This document is a chronological experiment log,
> and several conclusions further down have been overturned by later
> experiments. This section and [`README.md`](README.md) take precedence.

## Conclusions and Preconditions (front matter)

### The proposal

**Multi-axis grid partitioning.** At one node, run a single 1-D SAH sweep per
coordinate axis to obtain one global threshold per axis, then make a single
linear pass turning the per-axis comparisons into a bit key that names the
child. This yields a 2×2 (`Ours4`) or 2×2×2 (`Ours8`) grid partition.
**No binary tree is built and no collapse pass is run** — the wide node is
committed to directly. This is the structural difference from the existing
wide BVHs of Dammertz / Wald / Ylitie, all of which derive the wide tree
*after* a complete binary construction.

**The hybrids** recurse with an inline-widened binned top down to `N ≤ T`,
then switch to `Ours` beneath. The top performs `R` successive binned rounds
without recording the intermediate nodes, giving a `2^R`-ary node: `R = 2`
for 4-ary, `R = 3` for 8-ary.

### The claim (at **matched arity** and a **stated SIMD width**)

The hybrids (`Binned8+Ours8`, `Binned4+Ours4`) build **faster** than the
corresponding collapse baseline, with SAH **better in 2D (4.6–10 %)** and
**roughly equal in 3D (−0.8 % to +2.2 %)**.

Build-time evidence, stated by strength:

- **3D: solid.** −15.2 % on both meshes. Dragon 630.1 ± 7.5 vs 742.8 ± 7.3 ms
  (10.7 SD apart); Bunny 92.2 ± 1.9 vs 108.7 ± 1.0 ms (7.7 SD).
- **2D: directionally consistent but noisy.** All four configurations favour
  the proposal (−14.8 % to −23.0 %), but two separate by only ~1.5–1.7 SD.
  Run-to-run variation between invocations of the same binary reached 7 %,
  larger than some of the effects being claimed. **The 2D build-time advantage
  should be treated as suggestive, not established.**

SAH values are deterministic and reproduce bit-for-bit across runs; the
uncertainty above concerns timing only.

### Scope: the proposal only pays off at the bottom of a hybrid

`Ours4` / `Ours8` are **not** competitive as standalone builders, and this is
not a tuning problem — it is structural.

Per node, `Ours` does the same work as an exact sweep: three sorts and
`3(n-1)` candidate evaluations. Binned SAH does `O(n + B)` with 45
evaluations and **no sort at all**. So `Ours` reduces the *number* of split
levels (Dragon: max depth 25 -> 15 for `Ours4`) while raising the cost of
*each* level by more than the depth saving returns:

| Method | Per-node work | Evaluations | Max depth | Build [ms] |
|---|---|---|---|---|
| 2-way Sweep | 3 sorts, `O(n log n)` | `3(n-1)` | 24 | 2115 |
| Binned (B=16) | binning, `O(n+B)`, no sort | 45 | 25 | **682** |
| Ours4 | 3 sorts, `O(n log n)` | `3(n-1)` | **15** | 1370 |
| Ours8 | 3 sorts, `O(n log n)` | `3(n-1)` | **11** | 1029 |

`Ours4` cuts depth by 40 % and is still **2x slower than binned**. Counting
split rounds without pricing each round is exactly the kind of cost
evaluation this repository got wrong before.

The depth argument holds only where the sort is cheap — at the bottom of a
hybrid, where `n <= T = 32`. That is where the measured gain comes from:
`Binned8 inline` (binned all the way down) 722.7 ms -> `Binned8+Ours8` (Ours
below T) **630.1 ms**, a 13 % saving. Applied to a whole tree the same
argument fails.

**So the claim is scoped to: `Ours` as the bottom stage of a hybrid whose top
is an inline-widened binned build of matching arity.** Nothing broader is
claimed, and the standalone numbers above are not a fallback position.

### Retracted claim

An earlier version of this repository claimed a **"40–60 % SAH reduction."**
**That claim is wrong and is hereby withdrawn.** It came from using a model
that charges one traversal cost per inner node regardless of arity to compare
an 8-ary proposal against a 4-ary baseline. Nearly all of that gain was the
effect of *widening the tree*, which the existing collapse method (Dammertz
et al.) also achieves — and achieves slightly better. Every passage below
still asserting 40–60 % is superseded.

### What is explicitly not shown

- **No real ray traversal time was measured.** All SAH figures are proxy
  metrics. The build-time advantage is measured wall-clock time; the SAH
  comparison is a model.
- **The standalone proposals (`Ours4` / `Ours8`) are strictly dominated** by
  the collapse baseline — worse SAH *and* slower to build. Only the hybrids
  hold up.
- **The finished tree is larger**: 26 % more nodes in 2D, 47 % in 3D. Peak
  memory *during* construction is roughly 40 % lower (no intermediate binary
  tree is materialised), but the resulting tree costs more to store.
- **2D and 3D disagree on the SAH direction** (2D consistently better, 3D
  roughly level). This is currently unexplained.

---

## Preconditions

Every number here depends on the following assumptions. They are stated up
front because the central error in the earlier version of this work was
leaving assumption 1 unstated.

### 1. SIMD width — the load-bearing assumption

The whole-tree SAH cost is

```
SAH = C_trav · Σ_{n ∈ Inner} M(n)/M(root)  +  C_isect · Σ_{n ∈ Leaf} N(n)·M(n)/M(root)
```

where `M` is the half-perimeter (2D) or half surface area (3D) of a node's
box and `N(n)` is the primitive count of leaf `n`.

This form charges **one** `C_trav` per inner node **regardless of its arity
`k`**. That is correct only when the SIMD width `W` is at least `k`, so that
testing `k` child boxes is a single vector operation. It is the standard
model in the wide-BVH literature (Dammertz et al. 2008; Ernst & Greiner
2008) — but it is a model **with a precondition**, not a neutral metric.

**Consequence: SAH values from trees of different arity are not comparable
unless `W` is stated.** `sah_cost_simd(nodes, W)` charges `ceil(k/W)` vector
box tests per `k`-ary node and is the honest comparison; setting
`W ≥ max arity` recovers the flat form above.

The benchmarks print `maxK` (the largest arity actually produced) so
`W ≥ k` can be **checked** rather than assumed. All results here satisfy it:
4-ary methods reach `maxK = 4` and are compared at `W = 4`; 8-ary methods
reach `maxK = 8` and are compared at `W = 8`.

**At the wrong width the ranking inverts.** On Dragon, `Binned8+Ours8`
scores 24.32 at `W = 8` but 36.07 at `W = 4`, where the 4-ary
`Binned Collapse k=2` (30.48) overtakes it. **An 8-ary method has no
advantage on 4-wide hardware.**

### 2. Fixed parameters

| Parameter | Value |
|---|---|
| `C_TRAV`, `C_ISECT` | 1.0, 1.0 |
| `MAX_LEAF` | 4 primitives |
| Bin count `B` | 16 |
| Hybrid threshold `T` | 32 (runtime-settable with `--T=<n>`) |
| Area proxy | half-perimeter (2D) / half surface area (3D) |

### 3. Timing method

Untimed warm-up builds run first, then `n` timed builds; mean and **sample
standard deviation** are reported. 3D mesh runs use `n = 7` (`--runs=<n>`),
2D uses `n = 10`. Differences smaller than a few SD are not treated as real.
Timings are single-threaded and machine-specific; only *relative* differences
within the same run are meaningful.

### 4. Arity must be matched before comparing

A hybrid with a 4-ary top and an 8-ary bottom is not a fair opponent for a
uniformly 8-ary tree — SAH weight concentrates near the root, so a narrow top
wastes lanes exactly where cost is highest. Changing the top from two rounds
to three (`Binned4+Ours8` → `Binned8+Ours8`) moved Dragon's SAH from 29.60 to
24.32. **Most of the previously reported gap was this mismatch, not the
method.**

---

## Experimental Setup

| Parameter | Value |
|---|---|
| Dimension | 2D |
| Area proxy | Half-perimeter |
| C_TRAV | 1.0 |
| C_ISECT | 1.0 |
| MAX_LEAF | 4 (max primitives per leaf) |
| HYBRID_THRESH | 32 (hybrid strategy switch threshold) |
| Timing | Average of 10 runs per strategy |

### Test Primitives

- **Circles**: radius R=1, centers uniformly random
- **Triangles**: centered at random points, vertex offsets ∈ [-R, R]² (small random triangles)

### Test Scales

| N | range | Note |
|---|---|---|
| 128 | ±12 | Small scale, all strategies run |
| 1024 | ±36 | Medium scale, C/D strategies run |
| 102400 | ±340 | Large scale, C/D skipped |

---

## Strategy List

### 0. Median Split

**Description**: Bisect primitives at the centroid-range midpoint of the longest axis. No SAH; purely spatially balanced splitting.

- Arity: 2
- Time complexity: O(N log N) (sort per level)
- Heavyweight: no

---

### 1. 2-way SAH Sweep

**Description**: Sort by centroid, then forward/backward prefix-AABB sweeps evaluate all N-1 split positions; pick the SAH-minimal one. Both X and Y axes evaluated, best taken.

- Arity: 2
- Time complexity: O(N log N) (sort-dominated)
- Heavyweight: no

---

### 2. Binned SAH (B=16)

**Description**: Divide the centroid range into B=16 equal-width bins; evaluate only B-1=15 bin boundaries as split candidates, with prefix/suffix AABB sweeps. About 2-3× faster than Sweep SAH with <0.5% SAH loss.

- Arity: 2
- Time complexity: O(N + B) per node
- Heavyweight: no
- Note: when N < B many bins are empty, but behavior remains correct (empty bins skipped)

---

### 3. Collapse k=2 (Fast Quad-BVH)

**Description**: Build a full 2-way SAH Sweep tree, then post-process: each internal node skips one level and adopts its grandchildren as direct children, discarding the middle layer. Eliminates the middle layer's C_TRAV; result is roughly a 4-ary tree.

- Arity: ~4
- Time complexity: O(N log N) (build) + O(N) (collapse)
- Heavyweight: no
- Equivalent published method: Fast Quad-BVH (Dammertz et al.)

---

### 4. Binned Collapse k=2

**Description**: Same as Collapse k=2, but the build phase uses Binned SAH (B=16) for speed. SAH nearly identical to Collapse k=2; faster build.

- Arity: ~4
- Time complexity: O(N·B) (build) + O(N) (collapse)
- Heavyweight: no

---

### 5. A Independent (independent cross-axis combination)

**Description**: Find the best X-axis split (sx) and the best Y-axis split (sy). Use the cheaper axis as the primary split, then apply the other axis independently to each side, producing up to 4 children (2×2 cross-axis). **Same-axis secondary splits are not allowed.**

- Arity: 2 – 4
- Time complexity: O(N log N) (3 one-dimensional sweeps per node)
- Heavyweight: no

---

### 6. B Hierarchical (hierarchical cross-axis combination)

**Description**: First cut uses best_any (best of X/Y); second cut applies best_any again to each side, up to 4 children. Unlike A, secondary splits may reuse the primary axis.

- Arity: 2 – 4
- Time complexity: O(N log N) (3 one-dimensional sweeps per node)
- Heavyweight: no

---

### 7. C Exhaustive 2D

**Description**: For all (N-1)² combinations of X-order index sx and Y-order index sy, grid-partition and evaluate every cross-axis 4-way split; pick the SAH-minimal (sx, sy).

- Arity: 1 – 4 (empty cells merged)
- Time complexity: O(N³)
- Heavyweight: **yes** (skipped for N>1024)

---

### 8. D All-Axis

**Description**: Extends strategy C with same-axis double-cut candidates (X into Left/Mid/Right; Y into Left/Mid/Right), using incremental middle-segment AABBs to avoid recomputation.

- Arity: 1 – 4
- Time complexity: O(N³) (cross-axis) + O(N²) (same-axis) ≈ O(N³)
- Heavyweight: **yes** (skipped for N>1024)
- Finding: on uniformly random data, same-axis double cuts are almost never better than cross-axis; SAH nearly identical to C

---

### 9. Binned+A (T=32) (hybrid: Binned top + A bottom)

**Description**: Hybrid. While node primitive count N > 32, split 2-way with Binned SAH (B=16) and recurse; when N ≤ 32, switch to strategy A (independent cross-axis). Fast approximation on top, precise multi-way at the bottom.

- Arity: 2 on top, 2–4 at bottom
- Heavyweight: no
- Threshold: T = 32

---

### 10. Binned+SAH (T=32) (hybrid: Binned top + Sweep bottom)

**Description**: As above, but the bottom switches to 2-way SAH Sweep (exact binary). Control group for isolating the bottom strategy's effect.

- Arity: 2 throughout
- Heavyweight: no
- Threshold: T = 32

---

### 11. Binned4+A (T=32) (hybrid: 4-wide Binned top + A bottom)

**Description**: On top of Binned+A, each top node performs **two rounds** of inline Binned splitting (no intermediate nodes recorded), directly yielding up to 4 children — equivalent to collapsing the top. Bottom remains strategy A.

- Arity: ~4 on top (2×2), 2–4 at bottom
- Heavyweight: no
- Threshold: T = 32 (second round also checks the threshold)

---

### 12. Binned+SAH+Coll (T=32) (hybrid: Binned top + Sweep bottom + full-tree collapse)

**Description**: Build the Binned+SAH hybrid tree (binary throughout), then apply collapse_k2 to the whole tree, making it ~4-ary everywhere. Equivalent to "Binned Collapse k=2 with Binned replacing the top-level Sweep."

- Arity: ~4 throughout
- Heavyweight: no
- Threshold: T = 32 (controls the Binned/Sweep switch point)

---

## Results

> N=102400 differentiates the strategies best and is the primary reference.

### Circle Primitives

#### N=128

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 0.015 | 36.36 | 63 | 32 | 5 |
| 2-way SAH Sweep | 0.053 | 31.09 | 83 | 42 | 6 |
| Binned SAH (B=16) | 0.046 | 31.07 | 83 | 42 | 7 |
| Binned Collapse k=2 | 0.059 | 24.65 | 62 | 42 | 4 |
| Collapse k=2 | 0.056 | 24.65 | 62 | 42 | 3 |
| A Independent | 0.039 | 22.93 | 71 | 51 | 3 |
| B Hierarchical | 0.048 | 22.86 | 71 | 51 | 3 |
| C Exhaustive 2D | 5.47 | 21.73 | 81 | 61 | 4 |
| D All-Axis | 5.40 | 21.73 | 81 | 61 | 4 |
| Binned+A (T=32) | 0.039 | 26.09 | 71 | 48 | 5 |
| Binned+SAH (T=32) | 0.035 | 31.07 | 83 | 42 | 7 |
| Binned4+A (T=32) | 0.045 | 24.63 | 69 | 48 | 4 |
| Binned+SAH+Coll (T=32) | 0.075 | 24.65 | 62 | 42 | 4 |

#### N=1024

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 0.271 | 109.83 | 511 | 256 | 8 |
| 2-way SAH Sweep | 0.679 | 95.36 | 687 | 344 | 11 |
| Binned SAH (B=16) | 0.493 | 95.54 | 693 | 347 | 10 |
| Binned Collapse k=2 | 0.554 | 73.25 | 520 | 347 | 5 |
| Collapse k=2 | 0.791 | 73.31 | 516 | 344 | 6 |
| A Independent | 0.631 | 69.36 | 575 | 402 | 6 |
| B Hierarchical | 0.796 | 69.24 | 573 | 401 | 6 |
| C Exhaustive 2D | 2680 | 65.82 | 679 | 509 | 5 |
| D All-Axis | 2705 | 65.81 | 680 | 509 | 5 |
| Binned+A (T=32) | 0.468 | 76.53 | 666 | 445 | 9 |
| Binned+SAH (T=32) | 0.350 | 95.38 | 699 | 350 | 10 |
| Binned4+A (T=32) | 0.424 | 69.53 | 644 | 445 | 6 |
| Binned+SAH+Coll (T=32) | 0.433 | 72.87 | 523 | 350 | 5 |

#### N=102400

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 72.5 | 1123.5 | 65535 | 32768 | 15 |
| 2-way SAH Sweep | 144.4 | 1013.5 | 68251 | 34126 | 18 |
| Binned SAH (B=16) | 61.8 | 1018.4 | 68385 | 34193 | 18 |
| Binned Collapse k=2 | 67.0 | 788.1 | 52538 | 34193 | 9 |
| Collapse k=2 | 148.5 | 785.5 | 52510 | 34126 | 9 |
| A Independent | 128.4 | 715.1 | 62009 | 43508 | 10 |
| B Hierarchical | 148.6 | 710.4 | 62269 | 43885 | 9 |
| C Exhaustive 2D | SKIPPED | — | — | — | — |
| D All-Axis | SKIPPED | — | — | — | — |
| Binned+A (T=32) | 55.2 | 823.5 | 64821 | 43284 | 16 |
| Binned+SAH (T=32) | 49.9 | 1018.2 | 68285 | 34143 | 18 |
| Binned4+A (T=32) | 59.4 | 724.6 | 62103 | 43284 | 10 |
| Binned+SAH+Coll (T=32) | 58.2 | 788.5 | 52505 | 34143 | 9 |

---

### Triangle Primitives

#### N=128

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 0.015 | 31.48 | 63 | 32 | 5 |
| 2-way SAH Sweep | 0.038 | 26.45 | 87 | 44 | 7 |
| Binned SAH (B=16) | 0.047 | 26.16 | 85 | 43 | 7 |
| Binned Collapse k=2 | 0.097 | 20.18 | 64 | 43 | 4 |
| Collapse k=2 | 0.046 | 20.36 | 65 | 44 | 4 |
| A Independent | 0.035 | 17.95 | 77 | 57 | 4 |
| B Hierarchical | 0.044 | 17.10 | 82 | 61 | 4 |
| C Exhaustive 2D | 5.55 | 17.06 | 84 | 63 | 4 |
| D All-Axis | 5.57 | 17.06 | 84 | 63 | 4 |
| Binned+A (T=32) | 0.074 | 18.88 | 88 | 64 | 5 |
| Binned+SAH (T=32) | 0.035 | 26.16 | 85 | 43 | 7 |
| Binned4+A (T=32) | 0.070 | 17.40 | 86 | 64 | 4 |
| Binned+SAH+Coll (T=32) | 0.042 | 20.18 | 64 | 43 | 4 |

#### N=1024

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 0.289 | 93.75 | 511 | 256 | 8 |
| 2-way SAH Sweep | 0.665 | 77.58 | 693 | 347 | 10 |
| Binned SAH (B=16) | 0.539 | 78.62 | 695 | 348 | 10 |
| Binned Collapse k=2 | 0.493 | 58.40 | 524 | 348 | 5 |
| Collapse k=2 | 0.748 | 57.63 | 521 | 347 | 5 |
| A Independent | 0.615 | 51.24 | 635 | 464 | 5 |
| B Hierarchical | 0.792 | 51.07 | 640 | 466 | 5 |
| C Exhaustive 2D | 2693 | 50.40 | 687 | 514 | 6 |
| D All-Axis | 2704 | 50.40 | 687 | 514 | 6 |
| Binned+A (T=32) | 0.408 | 59.98 | 695 | 485 | 8 |
| Binned+SAH (T=32) | 0.458 | 78.50 | 695 | 348 | 10 |
| Binned4+A (T=32) | 0.453 | 52.62 | 670 | 485 | 6 |
| Binned+SAH+Coll (T=32) | 0.449 | 58.17 | 522 | 348 | 5 |

#### N=102400

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 69.8 | 934.9 | 65535 | 32768 | 15 |
| 2-way SAH Sweep | 138.7 | 822.0 | 68703 | 34352 | 18 |
| Binned SAH (B=16) | 59.1 | 826.4 | 68785 | 34393 | 18 |
| Binned Collapse k=2 | 66.4 | 617.8 | 52833 | 34393 | 9 |
| Collapse k=2 | 145.7 | 615.2 | 52824 | 34352 | 9 |
| A Independent | 129.6 | 516.5 | 69049 | 50590 | 9 |
| B Hierarchical | 149.3 | 512.5 | 69451 | 50979 | 9 |
| C Exhaustive 2D | SKIPPED | — | — | — | — |
| D All-Axis | SKIPPED | — | — | — | — |
| Binned+A (T=32) | 54.2 | 619.7 | 72193 | 50668 | 16 |
| Binned+SAH (T=32) | 51.4 | 826.2 | 68687 | 34344 | 18 |
| Binned4+A (T=32) | 58.0 | 524.4 | 69475 | 50668 | 10 |
| Binned+SAH+Coll (T=32) | 58.3 | 618.1 | 52803 | 34344 | 9 |

---

## Key Conclusions (N=102400, triangles)

### SAH ranking (worse → better)

```
Median > 2way SAH ≈ Binned SAH > Binned+SAH ≈ Binned+SAH+Coll ≈ Binned Collapse
  > Binned4+A > A Independent ≈ B Hierarchical > [C/D exhaustive, skipped]
```

### Speed vs quality trade-off

| Scheme | Time (ms) | SAH (tri, N=102400) | Reference |
|---|---|---|---|
| Binned SAH | 59 | 826 | Fast baseline |
| Binned Collapse k=2 | 66 | 618 | −25% SAH, +12% time |
| **Binned+SAH+Coll** | 58 | 618 | −25% SAH, ≈same speed |
| A Independent | 130 | 516 | Best non-exhaustive SAH, 2× slower |
| **Binned4+A** | 58 | 524 | Within 1.5% of A's SAH, 2× faster |
| B Hierarchical | 149 | 513 | Best non-exhaustive SAH, slowest |

### Takeaways

1. **Binned+SAH+Coll** is practically equivalent to Binned Collapse k=2; interchangeable implementations.
2. **Binned4+A** is the best practical scheme in this experiment: nearly A-Independent quality (524 vs 516, 1.5% gap) at roughly Binned-SAH build cost, 2× faster than A.
3. D All-Axis shows no meaningful gain over C Exhaustive on uniform random data — same-axis double cuts rarely help this distribution.
4. Differences shrink at small N (128); collapse's SAH benefit grows with scale.

---

## 3D Extension Experiment

### Setup

| Parameter | Value |
|---|---|
| Dimension | 3D |
| SAH proxy | Half surface area = dx·dy + dy·dz + dz·dx |
| C_TRAV | 1.0 |
| C_ISECT | 1.0 |
| MAX_LEAF | 4 |
| HYBRID_THRESH | 32 |
| Timing (synthetic) | Average of 10 runs |
| Timing (meshes) | N<200K: 5 runs, N<600K: 3 runs, N≥600K: 1 run |

### Test Primitives

- **Spheres**: radius R=1, centers uniformly random
- **Triangles (3D)**: centered at random points, vertex offsets ∈ [-R, R]³
- **Real meshes**: Stanford Bunny (144K triangles), Chinese Dragon (871K triangles)

---

### 3D Strategy List

#### Ported strategies (from 2D; `best_any` now enumerates 3 axes)

| Strategy | Description |
|---|---|
| Median Split | Longest-axis centroid-midpoint bisection |
| 2-way SAH Sweep | Sweep all 3 axes, take SAH-minimal split |
| Binned SAH (B=16) | 3-axis binning, B-1=15 candidates |
| Collapse k=2 | Sweep tree + collapse, ~4-ary |
| Binned Collapse k=2 | Binned tree + collapse, ~4-ary |

#### New 3D strategies

**A4 Independent (4-way independent cross-axis)**

- Primary: pick the SAH-cheapest of X/Y/Z for the first cut → L, R
- Secondary: L and R each independently pick the best of the **remaining 2 axes**
- Up to 4 children (2 per half)
- Arity: 2 – 4

**A8 Independent (8-way globally consistent axis order)**

- Primary: globally best of 3 axes → L, R
- Secondary: one axis chosen for both L and R (minimizing cost_L + cost_R) → up to 4 groups
- Tertiary: the remaining third axis applied to each group → up to 8 children
- Arity: 2 – 8; axis order consistent within a node (primary → secondary → tertiary)

**B4 Hierarchical (hierarchical 3-axis 4-way)**

- 3 greedy best_any cuts: first → L/R; second L → LL/LR, R → RL/RR
- Secondary cuts may reuse the primary axis
- Arity: 2 – 4

**B8 Hierarchical (hierarchical 3-axis 8-way)**

- 4th-level greedy: each of B4's four groups gets one more best_any cut
- Arity: 2 – 8

**Ours4 (global two-axis grid)**

- Run 3 one-dimensional SAH sweeps on the whole set; pick the 2 cheapest axes
- Use those 2 **global** split thresholds to form a 2×2 grid; assign primitives by centroid
- Key difference vs A4: the secondary threshold is shared by both halves (not adaptive)
- Arity: 2 – 4

**Ours8 (global three-axis grid)**

- Same, using all three axes' global optimal thresholds → 2×2×2 = 8 cells
- Arity: 2 – 8

#### Hybrid strategies (3D)

| Strategy | Description |
|---|---|
| Binned+A4 (T=32) | Binned 2-way top; A4 bottom (N≤32) |
| Binned4+A4 (T=32) | 2-round inline Binned top (~4-way); A4 bottom |
| Binned4+A8 (T=32) | 2-round inline Binned top (~4-way); A8 bottom |
| Binned4+Ours4 (T=32) | 2-round inline Binned top (~4-way); Ours4 bottom |
| Binned4+Ours8 (T=32) | 2-round inline Binned top (~4-way); Ours8 bottom |

---

### Synthetic Results (3D, N=102400)

#### Sphere primitives

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 77.3 | 155.20 | 65535 | 32768 | 15 |
| 2-way SAH Sweep | 224.2 | 133.05 | 67823 | 33912 | 17 |
| Binned SAH (B=16) | 77.5 | 134.10 | 67815 | 33908 | 17 |
| Collapse k=2 | 229.8 | 82.51 | 52706 | 33912 | 9 |
| Binned Collapse k=2 | 81.1 | 83.05 | 52514 | 33908 | 9 |
| A4 Independent | 220.5 | 67.27 | 74680 | 55987 | 9 |
| **A8 Independent** | 194.7 | **54.07** | 62028 | 52775 | 6 |
| B4 Hierarchical | 242.8 | 66.83 | 75079 | 56285 | 9 |
| **B8 Hierarchical** | 239.6 | **52.88** | 59900 | 51020 | 6 |
| Ours4 (best 2-axis) | 165.1 | 75.92 | 65856 | 46603 | 10 |
| Ours8 (best 3-axis) | 120.2 | 55.49 | 64637 | 52997 | 7 |
| Binned+A4 (T=32) | 81.8 | 101.97 | 77925 | 56228 | 15 |
| Binned4+A4 (T=32) | 82.8 | 68.78 | 75202 | 56228 | 9 |
| **Binned4+A8 (T=32)** | **80.5** | **66.45** | 61426 | 51561 | 9 |
| Binned4+Ours4 (T=32) | 69.2 | 76.56 | 66049 | 46533 | 10 |
| **Binned4+Ours8 (T=32)** | **61.6** | 67.76 | 64299 | 52060 | 10 |

#### Triangle (3D) primitives

| Strategy | Time (ms) | SAH | Nodes | Leaves | Max Depth |
|---|---|---|---|---|---|
| Median Split | 74.0 | 145.11 | 65535 | 32768 | 15 |
| 2-way SAH Sweep | 218.5 | 124.15 | 67749 | 33875 | 17 |
| Binned SAH (B=16) | 75.8 | 125.10 | 67619 | 33810 | 17 |
| Collapse k=2 | 227.2 | 75.59 | 52614 | 33875 | 9 |
| Binned Collapse k=2 | 85.3 | 76.12 | 52429 | 33810 | 9 |
| A4 Independent | 234.3 | 61.94 | 74924 | 56174 | 9 |
| **A8 Independent** | 191.9 | **49.15** | 62303 | 53024 | 6 |
| B4 Hierarchical | 237.6 | 61.79 | 74883 | 56144 | 9 |
| **B8 Hierarchical** | 236.5 | **48.05** | 60169 | 51263 | 6 |
| Ours4 (best 2-axis) | 166.7 | 69.29 | 66269 | 47016 | 10 |
| Ours8 (best 3-axis) | 119.8 | 50.81 | 64260 | 52737 | 7 |
| Binned+A4 (T=32) | 78.7 | 96.13 | 78076 | 56337 | 15 |
| Binned4+A4 (T=32) | 80.5 | 63.64 | 75352 | 56337 | 9 |
| **Binned4+A8 (T=32)** | **77.5** | **61.49** | 60882 | 51145 | 9 |
| Binned4+Ours4 (T=32) | 74.4 | 70.27 | 66407 | 46890 | 10 |
| **Binned4+Ours8 (T=32)** | **67.6** | 62.77 | 63945 | 51884 | 9 |

---

### Real-Mesh Results

#### Stanford Bunny (144,046 triangles)

| Strategy | Time (ms) | SAH |
|---|---|---|
| Median Split | 108.7 | 46.01 |
| 2-way SAH Sweep | 319.6 | 37.46 |
| Binned SAH (B=16) | 112.2 | 37.76 |
| Collapse k=2 | 325.2 | 23.28 |
| Binned Collapse k=2 | 118.2 | 23.44 |
| A4 Independent | 304.4 | 21.88 |
| A8 Independent | 272.1 | 17.67 |
| B4 Hierarchical | 326.5 | 21.55 |
| **B8 Hierarchical** | 354.9 | **16.49** |
| Ours4 (best 2-axis) | 219.4 | 24.63 |
| Ours8 (best 3-axis) | 160.5 | 19.22 |
| Binned+A4 (T=32) | 120.0 | 32.88 |
| Binned4+A4 (T=32) | 128.3 | 22.18 |
| **Binned4+A8 (T=32)** | **121.6** | **21.60** |
| Binned4+Ours4 (T=32) | 106.1 | 24.00 |
| **Binned4+Ours8 (T=32)** | **98.2** | 22.09 |

#### Chinese Dragon (871,306 triangles)

| Strategy | Time (ms) | SAH |
|---|---|---|
| Median Split | 747.0 | 67.20 |
| 2-way SAH Sweep | 2295.0 | 49.66 |
| Binned SAH (B=16) | 760.4 | 50.16 |
| Collapse k=2 | 2331.5 | 30.23 |
| Binned Collapse k=2 | 820.1 | 30.48 |
| A4 Independent | 2161.2 | 29.13 |
| A8 Independent | 1934.6 | 23.94 |
| B4 Hierarchical | 2352.9 | 28.75 |
| **B8 Hierarchical** | 2361.8 | **21.81** |
| Ours4 (best 2-axis) | 1476.6 | 32.38 |
| Ours8 (best 3-axis) | 1100.1 | 26.63 |
| Binned+A4 (T=32) | 764.7 | 45.31 |
| Binned4+A4 (T=32) | 814.1 | 29.86 |
| **Binned4+A8 (T=32)** | **789.2** | **28.64** |
| Binned4+Ours4 (T=32) | 706.8 | 31.52 |
| **Binned4+Ours8 (T=32)** | **688.5** | 29.60 |

---

### 3D Key Conclusions

> **⚠ Item 6 of this section is superseded.** Its claim that "Binned4+Ours8 is
> the fastest-building of all strategies" has been overtaken by `Binned8+Ours8`
> (Dragon, 630 ms), and its SAH comparison was made without stating a SIMD
> width and between trees of mismatched arity. See the front matter of this
> document and "Re-measurement at Matched Arity" below. Items 1–5 keep their
> relative ranking but likewise hold only when `W >= maxK`.


1. **A8 / B8 achieve the best SAH**: 8-way full-axis combination lowers SAH ~20% below the 4-way variants (Dragon: A8=23.9 vs A4=29.1). The third dimension provides split freedom 2D lacks.

2. **Binned4+A8 is the best engineering trade-off**: build time nearly matches Binned SAH (Dragon: 789 ms vs 760 ms), while SAH drops ~43% (28.6 vs 50.2), within 20% of pure A8's quality.

3. **2D conclusions hold in 3D**: Binned SAH ≈ Sweep SAH quality at 3× the speed; Collapse k=2 yields a large SAH drop (Dragon: ~40%); exhaustive methods (C/D) are infeasible in 3D at O(N³) and were not implemented.

4. **Real meshes vs random synthetic**: identical rankings — strategy ordering is robust to mesh type.

5. **Ours4 / Ours8 vs A4 / A8 (global grid vs adaptive splitting)**:
   - Ours uses global optimal 1D thresholds; secondary cuts are not adapted per subgroup. A4/A8 pick axes independently per group.
   - On Dragon, Ours4 SAH = 32.4 vs A4 = 29.1 — an ~**11%** gap; Ours8 = 26.6 vs A8 = 23.9, likewise ~**11%**.
   - Ours is ~30% faster than same-arity A/B (Dragon: Ours4 1477 ms vs A4 2216 ms) since no per-group axis re-sweeps are needed.
   - Ours8 is faster than Ours4 (1100 ms vs 1477 ms): three-axis cells hold fewer primitives, so recursion is shallower.
   - Conclusion: global grid splitting is a fast approximation of adaptive splitting with a fixed ~11% SAH penalty — suitable when build speed matters more than quality.

6. **Binned4+Ours hybrids (4-wide Binned top + Ours bottom)**:
   - **Binned4+Ours8 is the fastest builder of all strategies** (Dragon: 688 ms — ~10% faster even than pure Binned SAH's 760 ms) because the Ours8 bottom needs no per-group re-sweeps.
   - SAH is only ~3% worse than Binned4+A8 (Dragon: 29.6 vs 28.6) while ~13% faster.
   - Bunny agrees: 98 ms (fastest overall), SAH=22.1, within ~2% of Binned4+A8's 21.6.
   - Conclusion: for maximum build speed choose Binned4+Ours8 over Binned4+A8; the SAH gap (2–3%) is negligible.

---

## Follow-up Experiment: Warm-up Timing and a Sweep of the Threshold T

> The numbers in this section were re-measured **with warm-up** and are therefore not
> directly comparable to the tables above (which had none).

### Correction to the Timing Method

All earlier timings averaged in the first (cold) build, which made small-N results scatter by
up to ±30% — e.g. Binned+SAH on 2D circles at N=128 once measured 0.1008 ms and re-measured
at 0.035 ms. Both programs now perform **untimed warm-up builds** before timing: 2 for
synthetic scenes, 1 for the large meshes. After the fix, T-independent methods scatter by
only about **±1.5%** across 4 independent runs (Sweep on Dragon: 2168 / 2112 / 2117 / 2158 ms).

### The Threshold T Is Now a Runtime Parameter

`HYBRID_THRESH` is no longer a compile-time constant; it is set on the command line with
`--T=<n>`, and the printed strategy names show the actual value:

```
bvh_compute.exe   --T=256
bvh3d_compute.exe meshes/dragon/dragon.obj --T=256
```

### Newly Added Strategies

| Strategy | Dim | Description |
|---|---|---|
| Ours4 | 2D | One 1-D SAH sweep per axis, then a single 2×2 grid partition into 4 children |
| Binned+Ours4 (T) | 2D / 3D | Plain binned 2-way top; switches to Ours4 when N ≤ T (no collapse) |
| Binned4+Ours4 (T) | 2D / 3D | Two inline binned rounds per top node (≈4-ary); Ours4 at the bottom |
| **Binned4+Ours8 (T)** | 3D | Same top (≈4-ary); **Ours8** at the bottom: combines the best split of all three axes into a single 2×2×2 grid partition into 8 children |
| Binned+SAH (T) | 3D | Binned 2-way top; exact Sweep when N ≤ T (control group) |

> Ours8 needs three coordinate axes, so it exists only in 3D; the 2D tables have no Ours8 rows.

### Re-measured Results (warmed up, T=32)

#### 2D circles, N=102,400

| Strategy | Time (ms) | SAH |
|---|---|---|
| 2-way SAH Sweep | 138.4 | 1013.46 |
| Binned SAH (B=16) | 55.75 | 1018.37 |
| Ours4 (proposed 4-way) | 112.2 | 722.80 |
| Binned+Ours4 (T=32) | 47.83 | 825.95 |
| Binned Collapse k=2 (Fast Quad BVH) | 62.81 | 788.10 |
| Binned4+Ours4 (T=32) | 51.37 | 727.06 |
| Binned+SAH (T=32) (control) | 48.07 | 1018.20 |

#### 3D spheres, N=102,400

| Strategy | Time (ms) | SAH |
|---|---|---|
| 2-way SAH Sweep | 224.2 | 133.05 |
| Binned SAH (B=16) | 77.16 | 134.10 |
| Ours4 (proposed 4-way) | 158.6 | 75.92 |
| Ours8 (proposed 8-way) | 114.0 | 55.49 |
| Binned+Ours4 (T=32) | 65.97 | 109.75 |
| Binned Collapse k=2 (Fast Quad BVH) | 80.85 | 83.05 |
| Binned4+Ours4 (T=32) | 70.75 | 76.56 |
| Binned4+Ours8 (T=32) | 63.31 | 67.76 |
| Binned+SAH (T=32) (control) | 67.82 | 134.01 |

#### Stanford Bunny (144,046 triangles)

| Strategy | Time (ms) | SAH |
|---|---|---|
| 2-way SAH Sweep | 300.9 | 37.45 |
| Binned SAH (B=16) | 102.5 | 37.76 |
| Ours4 (proposed 4-way) | 207.9 | 24.63 |
| Ours8 (proposed 8-way) | 156.4 | 19.22 |
| Binned+Ours4 (T=32) | 98.29 | 34.70 |
| Binned Collapse k=2 (Fast Quad BVH) | 119.8 | 23.44 |
| Binned4+Ours4 (T=32) | 105.0 | 24.00 |
| Binned4+Ours8 (T=32) | 97.68 | 22.09 |
| Binned+SAH (T=32) (control) | 97.88 | 37.74 |

#### Chinese Dragon (871,306 triangles)

| Strategy | Time (ms) | SAH |
|---|---|---|
| 2-way SAH Sweep | 2168.5 | 49.66 |
| Binned SAH (B=16) | 695.0 | 50.16 |
| Ours4 (proposed 4-way) | 1445.4 | 32.38 |
| Ours8 (proposed 8-way) | 1096.1 | 26.63 |
| Binned+Ours4 (T=32) | 709.8 | 46.96 |
| Binned Collapse k=2 (Fast Quad BVH) | 773.1 | 30.48 |
| Binned4+Ours4 (T=32) | 803.8 | 31.52 |
| Binned4+Ours8 (T=32) | 721.8 | 29.60 |
| Binned+SAH (T=32) (control) | 717.5 | 50.16 |

### T Sweep (time ms / SAH)

#### 2D circles, N=102,400

| Strategy | T=32 | T=64 | T=256 | T=1024 |
|---|---|---|---|---|
| Binned+Ours4 | 47.83 / 825.95 | 49.44 / 814.24 | 53.63 / 782.78 | 63.39 / 765.62 |
| Binned4+Ours4 | 51.37 / 727.06 | 51.49 / 764.77 | 55.19 / 759.88 | 64.63 / 755.03 |
| Binned+SAH (control) | 48.07 / 1018.20 | 50.02 / 1018.22 | 57.40 / 1016.32 | 68.28 / 1014.55 |

#### 3D spheres, N=102,400

| Strategy | T=32 | T=64 | T=256 | T=1024 |
|---|---|---|---|---|
| Binned+Ours4 | 65.97 / 109.75 | 65.81 / 104.53 | 76.06 / 93.68 | 86.95 / 86.69 |
| Binned4+Ours4 | 70.75 / 76.56 | 68.67 / 83.89 | 78.33 / 81.53 | 91.82 / 79.98 |
| Binned4+Ours8 | 63.31 / 67.76 | 64.87 / 67.55 | 64.54 / 64.85 | 72.68 / 58.99 |
| Binned+SAH (control) | 67.82 / 134.01 | 71.73 / 134.02 | 88.46 / 133.79 | 106.1 / 133.40 |

#### Stanford Bunny

| Strategy | T=32 | T=64 | T=256 | T=1024 |
|---|---|---|---|---|
| Binned+Ours4 | 98.29 / 34.70 | 101.5 / 33.69 | 119.5 / 31.70 | 149.1 / 29.78 |
| Binned4+Ours4 | 105.0 / 24.00 | 108.6 / 24.01 | 117.3 / 24.23 | 132.5 / 24.35 |
| Binned4+Ours8 | 97.68 / 22.09 | 94.56 / 21.83 | 104.9 / 21.39 | 111.6 / 20.95 |
| Binned+SAH (control) | 97.88 / 37.74 | 101.0 / 37.73 | 125.2 / 37.66 | 156.9 / 37.55 |

#### Chinese Dragon

| Strategy | T=32 | T=64 | T=256 | T=1024 |
|---|---|---|---|---|
| Binned+Ours4 | 709.8 / 46.96 | 676.0 / 45.76 | 733.8 / 43.31 | 822.1 / 40.89 |
| Binned4+Ours4 | 803.8 / 31.52 | 712.9 / 31.71 | 756.8 / 32.07 | 880.3 / 32.35 |
| Binned4+Ours8 | 721.8 / 29.60 | 638.4 / 29.31 | 681.0 / 28.94 | 758.0 / 28.55 |
| Binned+SAH (control) | 717.5 / 50.16 | 677.2 / 50.15 | 782.4 / 50.11 | 946.5 / 50.01 |

### Conclusions from the T Sweep

> **⚠ Item 6 of this section is wrong and is withdrawn.** It claimed that
> "Binned4+Ours8 beats Fast Quad BVH on both time and SAH" (Dragon, T=64:
> 638 ms / 29.31 versus 735 ms / 30.48). That compared an **8-ary** proposal
> against a **4-ary** Fast Quad under a cost model insensitive to arity. At
> the matched width `W = 8`, Fast Quad's 8-ary counterpart
> (`Binned Collapse 8ary`) scores **23.81**, clearly better than 29.31. The
> remaining items — the monotonicity in T and the mechanism in item 7 — are
> unaffected, because they compare T within a **single** arity.


1. **Binned+Ours4: T is a clean quality/speed knob.** Raising T lowers SAH and raises build
   time, both monotonically (Dragon: 47.0 → 40.9, −13%; 710 → 822 ms, +16%). Even at T=1024,
   however, its SAH does not reach what Binned4+Ours4 already achieves at T=32.
2. **Binned4+Ours4: T=32 is already optimal; larger T is slightly worse** (Dragon: 31.5 → 32.4).
   The SAH gain comes mainly from merging intermediate nodes at the top, and a larger T means
   less of that merging — so there is no reason to tune T for this method.
3. **Binned+SAH (control): SAH is essentially insensitive to T** (Dragon: 50.16 → 50.01, a
   0.3% spread across four values of T). This confirms that switching the bottom to a more
   precise **binary** method yields nothing; the gain comes only from switching it to a
   **4-way** split.
4. **Binned4+Ours8 (three-axis combination, 3D only) beats Binned4+Ours4 on every 3D scene at
   every T — in both time and SAH.** Dragon at T=64: 638 ms / 29.31 against 712 ms / 31.71 for
   the Ours4 version; Bunny at T=64: 94.6 ms / 21.83 against 108.6 ms / 24.01. Splitting all
   three axes at once leaves fewer primitives per cell and a shallower tree (max depth 13 vs 15
   on Dragon), while avoiding the repeated recursion the Ours4 version performs at the bottom.

5. **For Binned4+Ours8 a larger T is better — the opposite of the Ours4 version.**
   Dragon's SAH falls monotonically with T: 29.60 → 29.31 → 28.94 → 28.55; Bunny: 22.09 → 21.83
   → 21.39 → 20.95; 3D spheres: 67.76 → 67.55 → 64.85 → 58.99. Binned4+Ours4, in contrast, gets
   worse beyond T=32. In other words, the 8-way split is worth handing more of the tree to; the
   4-way split is not.

6. **Binned4+Ours8 beats Fast Quad BVH (the established method) on time and SAH at once.**
   Dragon at T=64: 638 ms / 29.31 against 735 ms / 30.48;
   Bunny at T=64: 94.6 ms / 21.83 against 108.7 ms / 23.44;
   3D spheres at T=32: 63.3 ms / 67.76 against 79.3 ms / 83.05.
   It is so far the only scheme that dominates the established method on both metrics.

7. **The opposite responses of 4 and 6 to T have one explanation.** Raising T shrinks the
   binned top and enlarges the part of the tree built by the proposed rule, so the hybrid's SAH
   moves monotonically toward that of the corresponding *pure* method. Which direction depends
   only on whether the pure method beats the binned top:
   - On Dragon, pure Ours4 = 32.38, **worse** than the binned top, so Binned4+Ours4 degrades
     from 31.52 (T=32) toward 32.35 (T=1024), approaching 32.38.
   - Pure Ours8 = 26.63, **better** than the binned top, so Binned4+Ours8 improves from 29.60
     toward 28.55, likewise approaching 26.63.

   Both curves are monotone and both converge on their pure method's value, exactly as the
   mechanism predicts. Practical reading: confine the four-way rule to small nodes (small T),
   and give the eight-way rule as much of the tree as the build budget allows (large T).

---

## Re-measurement at Matched Arity (current valid conclusions)

> This section supersedes every cross-arity SAH comparison above. See
> "Preconditions" in the front matter for the method.

### Controls added

| Name | Description |
|---|---|
| `Binned Collapse 8ary` | Binned binary tree collapsed by **2** levels → 8-ary. The existing baseline `Ours8` should actually be compared against |
| `Binned4 inline (pure)` | Pure 2-round inline binned, no Ours bottom. Control for R=2 |
| `Binned8 inline (pure)` | Pure 3-round inline binned, no Ours bottom. Control for R=3 |
| `Binned8+Ours8` | 3-round inline binned top (8-ary) + Ours8 bottom. **Top and bottom arity matched** |

### Verified: inline widening ≡ collapsing

`R` inline binned rounds produce a tree **topologically identical** to
building a binary binned tree and collapsing it by `R−1` levels. SAH, node
count, leaf count and max depth match **exactly** — verified for `R = 2` and
`R = 3`, in 2D and 3D, in every test configuration.

| Dragon | SAH(W=8) | inner | nodes | leaves | maxdepth |
|---|---|---|---|---|---|
| `Binned Collapse 8ary` | 23.8106 | 98,157 | 385,715 | 287,558 | 9 |
| `Binned8 inline (pure)` | 23.8106 | 98,157 | 385,715 | 287,558 | 9 |

| 2D N=102400 circles | SAH(W=4) | SIMD2 | nodes |
|---|---|---|---|
| `Binned Collapse k=2` | 788.1023 | 962.0600 | 52,538 |
| `Binned4 inline (pure)` | 788.1023 | 962.0600 | 52,538 |

This makes the ablation valid: `Binned8+Ours8` versus `Binned Collapse 8ary`
have an **identical top and differ only in the bottom**, so the measured
difference is attributable to Ours. It also identifies the source of the
build-time saving — the intermediate binary tree is never materialised.

### 3D real meshes, `W = 8`, `T = 32`, n = 7

Dragon (871,306 triangles):

| Method | SAH (W=8) | Build [ms] | Nodes | maxK |
|---|---|---|---|---|
| Binned SAH (B=16), binary baseline | 50.163 | 681.6 ± 6.8 | 575,115 | 2 |
| Binned Collapse k=2 (existing, 4-ary) | 30.477 | 748.1 ± 4.3 | 431,685 | 4 |
| **Binned Collapse 8ary (existing, 8-ary)** | **23.811** | 742.8 ± 7.3 | 385,715 | 8 |
| Binned8 inline (pure, same topology) | 23.811 | 722.7 ± 5.7 | 385,715 | 8 |
| **Binned8+Ours8 (proposed)** | 24.323 (+2.2 %) | **630.1 ± 7.5 (−15.2 %)** | 565,372 | 8 |
| Binned4+Ours8 (mismatched top, old) | 29.596 | 669.4 ± 6.7 | 571,179 | 8 |
| Ours8 standalone | 26.633 | 1028.7 ± 7.4 | 553,901 | 8 |
| Ours4 standalone | 32.382 (W=4) | 1370.5 | 542,343 | 4 |

Bunny (144,046 triangles):

| Method | SAH (W=8) | Build [ms] | Nodes | maxK |
|---|---|---|---|---|
| Binned SAH (B=16), binary baseline | 37.759 | 101.6 ± 0.9 | 91,773 | 2 |
| Binned Collapse 8ary (existing) | 18.714 | 108.7 ± 1.0 | 62,865 | 8 |
| **Binned8+Ours8 (proposed)** | **18.560 (−0.8 %)** | **92.2 ± 1.9 (−15.2 %)** | 95,074 | 8 |
| Ours8 standalone | 19.220 | 150.2 ± 3.1 | 91,101 | 8 |

### 2D synthetic, `W = 4`, `T = 32`, n = 10

`Binned4+Ours4` relative to `Binned Collapse k=2` (negative is better):

| Configuration | ΔSAH | ΔBuild time | Timing separation |
|---|---|---|---|
| N=1024, circles | −4.6 % | −18.2 % | ~1.5 SD (marginal) |
| N=102400, circles | −7.7 % | −17.2 % | ~3.3 SD |
| N=1024, triangles | −5.7 % | −23.0 % | ~13 SD |
| N=102400, triangles | −10.0 % | −14.8 % | ~1.7 SD (marginal) |

N=102400 circles in detail:

| Method | SAH (W=4) | SIMD2 | Build [ms] | nodes |
|---|---|---|---|---|
| Binned SAH (B=16) | 1018.3672 | 1018.3672 | 56.31 ± 2.02 | 68,385 |
| Binned Collapse k=2 | 788.1023 | 962.0600 | 62.37 ± 2.66 | 52,538 |
| Binned4 inline (pure) | 788.1023 | 962.0600 | 59.71 ± 2.83 | 52,538 |
| **Binned4+Ours4 (proposed)** | **727.0625** | 956.9994 | **51.62 ± 1.91** | 66,097 |

> **Note**: across two independent invocations of the same configuration,
> `Binned Collapse k=2` measured 58.19 ms and then 62.37 ms — a 7 % shift.
> **Timings are not comparable across runs; only within-run relative
> differences are meaningful.** SAH is unaffected (deterministic, bit-identical).

### Conclusions of this section

1. **Build time**: 3D is −15.2 % with solid evidence (Dragon 10.7 SD, Bunny
   7.7 SD). 2D is −14.8 % to −23.0 %, consistent in direction across all four
   configurations, but two are only ~1.5–1.7 SD apart and run-to-run drift
   reaches 7 %, so 2D should be read as **suggestive, not established**.
2. **SAH is consistently better in 2D (4.6–10 %) and roughly level in 3D**
   (Bunny −0.8 %, Dragon +2.2 %). The two disagree in direction; unexplained.
3. **The standalone proposals are dominated**: `Ours4` 32.38 / 1370 ms versus
   `Binned Collapse k=2` 30.48 / 748 ms — worse SAH *and* 83 % slower. Same
   for `Ours8`. **Only the hybrids stand up.**
4. **Cost**: the finished tree carries 26 % more nodes in 2D, 47 % in 3D. Peak
   memory during construction is roughly 40 % lower (collapse holds the binary
   tree, 575k nodes, and the new tree, 386k, simultaneously; inline builds
   565k directly).
5. **The proposal does not own the SAH frontier**: `B8 Hierarchical` reaches
   better SAH (Dragon 21.81, Bunny 16.49) but builds 3× slower. The proposal
   occupies a point on the build-time/SAH trade-off, not the quality optimum.
6. **All conclusions presuppose `W ≥ maxK`**, and SAH is a proxy — real
   traversal time has not been measured.

### Open problems

1. **Measure real ray traversal time.** A 2 % SAH difference should not be
   argued from the model alone.
2. **Explain the 2D/3D disagreement in SAH direction.**
3. **Make the split decision match the structure built.** Thresholds are
   chosen by per-axis *binary* SAH, but a 4-way or 8-way node is then built;
   the cost of that wide split is never evaluated, nor compared against the
   binary alternative. The criterion and the constructed structure do not
   correspond.
4. **Account for the larger finished tree (26–47 % more nodes)** in any memory
   argument.
5. **Reduce 2D timing noise** (more runs, pinned clocks) so the 2D build-time
   claim can be settled.
