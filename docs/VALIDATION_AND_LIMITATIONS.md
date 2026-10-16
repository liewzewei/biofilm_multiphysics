# Verification, Validation and Limitations

This document assesses how far the structural model described in [METHODOLOGY.md](METHODOLOGY.md) can be
trusted. Every number has one of three provenances:

| Tag | Meaning |
|---|---|
| **(F)** | Recomputed by [`notebooks/FINAL_biofilm_mechanics.ipynb`](../notebooks/FINAL_biofilm_mechanics.ipynb) on every run |
| **(P)** | Previously obtained in exploratory notebook 02, 03 or 04 and quoted here unchanged |
| **(D)** | Diagnostic run with the `src/` modules while writing the documentation (Stage 3, or Stage 4 after the strain fix in §3.5). It is not in any notebook |

---

## 1. Verification versus validation

| Activity | Question it answers | Status |
|---|---|---|
| **Code verification** | Does the code solve the stated equations correctly? | **Done** for geometry, kinematics, boundary conditions, rigid-body invariance, linear scaling, the cohesive failure index and grid conservation (§2, §3, §9). One postprocessing error in the principal strains was found and fixed (§3.5). |
| **Numerical convergence** | Are results insensitive to mesh resolution? | **Done** for bulk stress and strain indicators and for interface tractions (§4, §8). |
| **Experimental validation** | Do the equations and parameters describe real biofilms on real bellows? | **Not done.** No measurements of our engineered biofilm's mechanics or adhesion, of its thickness under any proposed thinning intervention, or of bellows kinematics have been compared with the model. |

Verification and convergence show that the numbers are *correct solutions of the model*. Only validation could
show that the *model describes reality*. Every conclusion in §11 is conditional on the assumptions in §10.

---

## 2. Geometry verification

| Check | Result | Criterion | Source |
|---|---|---|---|
| Bellows arc length vs. analytical $2l+R(\pi-\theta)$ | relative error $1.7\times10^{-6}$ | $<10^{-4}$ | (F) |
| Total length kept, $60^\circ\rightarrow60.5^\circ$ | relative change $1.4\times10^{-8}$ | $<10^{-5}$ | (F) |
| Total length kept, $60^\circ\rightarrow90^\circ$ (illustrative) | relative change $7.6\times10^{-7}$ | not asserted | (F) |
| Reference crevice independent of final angle | max difference $0$ | $<10^{-8}$ μm | (F) |
| Material-coordinate correspondence: length of resampled polyline, reference vs. deformed | relative change $7.9\times10^{-8}$ | $<10^{-5}$ | (F) |
| Substrate mesh nodes = reference bellows points | max difference $0$ | $<10^{-8}$ μm | (F) |
| Interface edge normals point into the biofilm | 499/499 edges, both conditions | 100 % | (F) |

**Mesh orientation and quality** for the production mesh ($N=500$, $n_y=15$, 14 970 triangles) (F):

| Condition | Min signed area | Min angle | Polygon valid | Exposed surface simple | Mesh vs. polygon area |
|---|---:|---:|:-:|:-:|---:|
| Unregulated | 0.1348 μm² | 20.60° | ✓ | ✓ | $4.3\times10^{-16}$ |
| Regulated | 0.0743 μm² | 10.72° | ✓ | ✓ | $6.5\times10^{-15}$ |

Uniform arc-length resampling was essential. The raw, non-uniformly sampled profile gave minimum angles of about
7° at $n_y=20$, falling to 2.4° at $n_y=60$ (P).

**Remaining self-intersection limitations.** The checks rule out inverted elements and a self-intersecting
exposed surface. They do not prove that interior elements never overlap, and the mesher does not enforce the
offset condition $h\kappa<1$. That condition holds here: the largest baseline thickness is $\max h=15.5\ \mu$m,
below $R_0=20\ \mu$m. In the fold-radius study at $R_0=15\ \mu$m, $h\kappa<1$ is not guaranteed by construction;
validity there rests on the post-mesh checks, which passed. The closest approach between non-neighbouring points
of the exposed surface was 4.6 μm in an earlier non-uniform-mesh test (P); it was not re-evaluated for the
production mesh.

---

## 3. FEM verification

### 3.1 Zero displacement

With $\mathbf u_b=\mathbf 0$ and no loads, the exact solution is $\mathbf u=\mathbf 0$ and $\boldsymbol\sigma=\mathbf 0$.
On the seed-42 unregulated production mesh, `solve_curved_biofilm` returns $\max\lvert\mathbf u\rvert=0$ and
$\max\sigma_e=0$ **(F)**; the max principal strain is also $0$ (D).

> **Provenance note.** In `04_bellows_deformation.ipynb` the zero-displacement test code sits inside a
> **Markdown** cell (cell 52), so it was never executed there. It was first run in Stage 3 (D) and is now an
> asserted check in Section 10 of the final notebook.

### 3.2 Rigid translation

Prescribing $\mathbf u_b=(0.5,\,0.3)\ \mu$m at every substrate node should translate the whole biofilm without
stress (F):

- displacement field error $9.9\times10^{-21}$ m
- $\max\sigma_e=2.3\times10^{-13}$ Pa
- max principal strain $4.6\times10^{-15}$ (P)

### 3.3 Infinitesimal rigid rotation

For $\mathbf u=\omega(-y,\,x)$ with $\omega=10^{-3}$ rad, the displacement gradient is antisymmetric, so
$\boldsymbol\varepsilon=\mathbf 0$ (F):

- field error $1.9\times10^{-21}$ m
- $\max\sigma_e=6.9\times10^{-14}$ Pa
- max principal strain $1.3\times10^{-15}$ (P)

This verifies invariance to *infinitesimal* rotation only. A finite rotation is not strain-free in linear
elasticity (§10.6).

### 3.4 Boundary-condition accuracy, solver convergence and linearity

| Check | Result | Source |
|---|---|---|
| FEM interface displacement vs. prescribed $\mathbf u_b$, both conditions | max error $5.6\times10^{-17}$ μm | (F) |
| Newton/direct solve, relative residual after one step | $2.3$–$2.6\times10^{-15}$ | (D) |
| Stress scales with $E$ ($10\times E\Rightarrow10\times\sigma$) | ratio $10.000000$ | (F) |
| Flat substrate: uniform-stretch BC error / sinusoidal BC error | $5.6\times10^{-17}$ / $8.9\times10^{-16}$ μm | (P) |
| Flat substrate: stretch vs. compression antisymmetry of $\mathbf u$, $\boldsymbol\sigma$ | error $0$ | (P) |
| Flat substrate: global force balance (element-stress integration along base) | 0.506 % | (P) |
| Test suite: flat-solver regression, mesh tests, strain-convention tests (`tests/`) | 8/8 pass | (D) |

### 3.5 Known inconsistency: shear-strain convention

**Status: fixed in Stage 4.** The heading is kept so existing links still resolve.

**Problem.** SfePy's `ev_cauchy_strain` returns $[\varepsilon_{xx},\varepsilon_{yy},\gamma_{xy}]$ with
$\gamma_{xy}=2\varepsilon_{xy}$. This is documented in SfePy's `CauchyStrainTerm` and was confirmed by a
homogeneous pure-shear test. Before Stage 4, `solve_curved_biofilm` and `solve_biofilm` inserted that component
directly as $\varepsilon_{xy}$ in

$$
\varepsilon_{1,2}=\tfrac{\varepsilon_{xx}+\varepsilon_{yy}}{2}\pm\sqrt{\left(\tfrac{\varepsilon_{xx}-\varepsilon_{yy}}{2}\right)^2+\varepsilon_{xy}^2},
$$

so `max_principal_strain`, `principal_tensile_strain` and `max_tensile_strain` were overestimated wherever shear
was present.

**Fix.** Both solvers now use `exy = 0.5 * strain_elements[:, 2]`. The returned `result["strain"]` array is
unchanged: it still holds SfePy's engineering shear strain in column 2.

**Not affected:** the stress output, which contains the true $\tau_{xy}$. Von Mises stress, principal tensile
stress, interface tractions, grid mapping and every stress-based result are therefore identical before and
after the fix. The final notebook asserts this.

**Tests** (`tests/test_strain_convention.py`):

1. SfePy's third strain component equals $2\varepsilon_{xy}$ under homogeneous pure shear.
2. Principal strains from `solve_curved_biofilm` are invariant when the geometry and the prescribed
   displacement are rotated by 30°. The pre-fix code failed this test, giving a 47 % difference.
3. `solve_biofilm`'s `max_principal_strain` equals the tensor eigenvalues computed from its own strain output.

**Effect on reported strain results:**

| Quantity | Before fix | After fix |
|---|---:|---:|
| Seed-42 max principal strain, unregulated / regulated | 0.5946 % / 0.2708 % | 0.4493 % / 0.1732 % (F) |
| Seed-42 max strain reduction | 54.46 % | 61.44 % (F) |
| 20-seed mean (median) max strain reduction | 57.15 % (57.18 %) | 60.37 % (60.84 %) (F) |
| 20-seed range of max strain reduction; pairs reduced | 52.15–63.59 %; 20/20 | 53.38–67.33 %; 20/20 (F) |
| Flat-substrate baseline max principal strain (notebook 02 case) | 0.1273 % | 0.1135 % (D) |

Notebook 04 has since been rerun with the corrected solvers, and its strain outputs match the values above.
Notebooks 02 and 03 were not rerun and still show pre-fix strain values. All strain values
quoted elsewhere in this document are post-fix, unless they are explicitly marked as pre-fix.

### 3.6 Cohesive failure-index verification

The threshold postprocessing of [METHODOLOGY.md §10.5](METHODOLOGY.md#105-picioreanu-inspired-cohesive-vulnerability-threshold-postprocessing)
is checked in the final notebook, Section 16 (F):

| Check | Result | Criterion |
|---|---|---|
| Synthetic three-element test ($\sigma_e=0.05,0.10,0.20$ Pa, $A=1,2,3$): $f_{\mathrm{over}}$ at $\sigma_c=0.15$, 10, 0.001 Pa | 0.5, 0, 1 | exact |
| Tie case: $\sigma_c=0.10$ Pa, element with $FI=1$ exactly | not counted ($f=0.5$) | strict inequality $FI>1$ |
| Representative pair at $\sigma_c=0.1$ Pa vs. notebook 04 (P) | max $FI$ 2.694 / 1.063; $f_{\mathrm{over}}$ 1.78 % / 0.11 % | identical to rounding |
| Phase-diagram entry vs. direct calculation | identical ($\Delta f=0.000970$ at the tested grid point) | $<10^{-12}$ |
| Sign convention: $\Delta f$ at $\beta=1$ over all baseline strengths | $\max\Delta f=0$ | $\le0$ |
| Fractions vanish above the peak stress | both zero from $\sigma_c=0.2736$ Pa (grid value); baseline max $\sigma_e=0.2694$ Pa | consistent |

The analysis operates on $\sigma_e$ only, so it is unaffected by the Stage 4 strain-convention fix (§3.5).

---

## 4. Mesh convergence (curved FEM)

The mesh was refined simultaneously along $s$ ($N$) and through the thickness ($n_y$), on the seed-42 profile,
with the geometry and $\mathbf u_b$ interpolated to each resolution (P). Relative change between successive
meshes:
$\delta_Q=\lvert Q_{\text{fine}}-Q_{\text{coarse}}\rvert/\lvert Q_{\text{fine}}\rvert$.

**Maximum von Mises stress (Pa)**

| $(N,n_y)$ | Elements | Unregulated | Regulated |
|---|---:|---:|---:|
| (250, 8) | 3 984 | 0.259431 | 0.105112 |
| **(500, 15)** (production) | **14 970** | **0.269396** | **0.106293** |
| (750, 23) | 34 454 | 0.275188 | 0.107499 |
| (1000, 30) | 59 940 | 0.275607 | 0.107352 |

**Relative change $\delta_Q$ between successive meshes**

| Metric | Condition | 1→2 | 2→3 | 3→4 |
|---|---|---:|---:|---:|
| Max σ_e | unregulated | 3.70 % | 2.10 % | 0.15 % |
| | regulated | 1.11 % | 1.12 % | 0.14 % |
| Mean σ_e | unregulated | 0.77 % | 0.05 % | 0.04 % |
| | regulated | 1.88 % | 0.14 % | 0.10 % |
| Max principal strain† | unregulated | 0.63 % | 0.82 % | 0.55 % |
| | regulated | 0.45 % | 0.58 % | 0.52 % |

**Convergence of the thinning effect (reduction, %)**

| Elements | Max σ_e | Mean σ_e | Max strain |
|---:|---:|---:|---:|
| 3 984 | 59.484 | 54.544 | 61.025 |
| **14 970** | **60.544** | **55.037** | **61.443** |
| 34 454 | 60.936 | 55.074 | 61.535 |
| 59 940 | 61.049 | 55.099 | 61.522 |

† The strain rows were recomputed after the Stage 4 fix (D), on the same meshes and inputs. The stress rows are
(P) and unaffected by the fix.

**Interpretation.**

- Stress quantities change by less than 0.3 % between the two finest meshes.
- Post-fix maximum strain fluctuates by 0.5–0.8 % between successive meshes rather than decreasing
  monotonically, which is typical of a single-element peak. Its reduction is stable to within 0.1 percentage
  points.
- On the production mesh, the unregulated maximum stress is about 2.3 % below the finest-mesh value.
- The *reduction* converges faster than the individual peaks, because discretization error is similar for
  both members of a pair.
- On the production mesh, the max-stress reduction is underestimated by about 0.5 percentage points; the
  mean-stress and max-strain reductions are within 0.1 points of the finest mesh.
- Unlike the flat-substrate model, whose fixed-base corner peak did not converge (§4.1), the curved model's
  global maximum converges. The reported maxima include the strip ends, where the clamped interface meets a
  free boundary.

### 4.1 Flat-substrate convergence (P, notebook 02)

| Metric | 1→2 | 2→3 | 3→4 |
|---|---:|---:|---:|
| Global max σ_e (includes fixed-base corners) | 17.23 % | 10.68 % | 7.79 % |
| Interior max σ_e (5 μm corner exclusion) | 0.52 % | 0.17 % | 0.14 % |
| Max displacement | 1.60 % | 0.52 % | 0.26 % |

The flat model's corner peak behaves like a boundary-condition singularity. Only interior quantities were used
for its conclusions.

---

## 5. Stochastic robustness

20 paired morphologies (seeds 0–19), $60^\circ\rightarrow60.5^\circ$, $E=64$ Pa, $\nu=0.3$, production mesh.
Recomputed (F). Stress rows are identical to (P) within 0.004 percentage points; the strain row reflects the
Stage 4 fix (pre-fix value: 57.15 %).

| Reduction after 50 % thinning | Mean | Median | Min | Max | Pairs reduced |
|---|---:|---:|---:|---:|:-:|
| Max von Mises stress | **59.57 %** | 60.13 % | 53.07 % | 65.51 % | 20/20 |
| Max principal strain | **60.37 %** | 60.84 % | 53.38 % | 67.33 % | 20/20 |
| Area-weighted mean von Mises stress | **52.22 %** | 53.53 % | 42.61 % | 57.48 % | 20/20 |

**Representative pair (seed 42) (F)**

| Metric | Unregulated | Regulated | Reduction |
|---|---:|---:|---:|
| Max von Mises stress | 0.269396 Pa | 0.106293 Pa | 60.54 % |
| Mean von Mises stress | 0.011936 Pa | 0.005367 Pa | 55.04 % |
| Max principal tensile stress | 0.327963 Pa | 0.123915 Pa | **62.22 %** |
| Mean principal tensile stress | 0.013812 Pa | 0.005790 Pa | 58.08 % |
| Max principal tensile strain | 0.4493 % | 0.1732 % | 61.44 % |
| Max displacement | 0.8004 μm | 0.7917 μm | 1.09 % |

The displacement hardly changes because it is dominated by the prescribed, near-rigid motion of the flanks.

**What the ensemble establishes.** Under fixed material, geometry and kinematics, the sign and approximate size
of the thinning effect do not depend on the particular random arrangement of surface features: every pair is
reduced, and the spread across pairs (about 10–15 percentage points) is much smaller than the effect.

**What it does not establish:**

- robustness to uncertainty in $E$, $\nu$, bellows geometry or the degree of thinning ($\alpha$ is fixed at 0.5);
- that any biological intervention actually halves thickness;
- any statistical statement about real ISS biofilms. The 20 samples come from a synthetic generator, not from
  imaged biofilms.

**Earlier flat-substrate ensemble (P, notebook 03).** For 20 pairs on a flat substrate under equal total
surface-traction force, the interior maximum von Mises stress fell by a mean of 23.96 % (median 24.07 %), with
20/20 pairs reduced; maximum displacement fell by 55.53 %. That loading case is different and its numbers
should not be mixed with the bellows results.

---

## 6. Material-parameter sensitivity

### 6.1 Young's modulus (P, notebook 04; spot-checked (F))

| $E$ (Pa) | 32 | 64 | 128 | 640 | 6400 |
|---|---:|---:|---:|---:|---:|
| Max tensile stress, unregulated (Pa) | 0.163982 | 0.327963 | 0.655927 | 3.279633 | 32.796329 |
| Max tensile stress, regulated (Pa) | 0.061957 | 0.123915 | 0.247830 | 1.239148 | 12.391475 |

At every $E$ the reductions are identical: max tensile stress 62.2169 %, mean tensile stress 58.0790 % (P),
and max tensile strain 61.443 % (D, post-fix; 54.4609 % pre-fix). Stress divided by $E$ was constant to $10^{-5}$ relative. This is exactly what
linearity predicts (METHODOLOGY §12.1). The **relative** conclusions therefore do not depend on the
uncalibrated $E$; **absolute** stresses do, in direct proportion.

### 6.2 Poisson's ratio (P; recomputed (F))

| $\nu$ | 0.20 | 0.30 | 0.40 | 0.45 |
|---|---:|---:|---:|---:|
| Max tensile stress, unregulated (Pa) | 0.2992 | 0.3280 | 0.3906 | 0.4752 |
| Max tensile stress reduction | 61.63 % | 62.22 % | 63.75 % | 66.00 % |
| Mean tensile stress reduction | 56.69 % | 58.08 % | 60.18 % | 61.73 % |
| Max tensile strain reduction (D, post-fix) | 61.30 % | 61.44 % | 61.74 % | 61.94 % |

Thinning reduces loading at every $\nu$ tested, and the reduction grows with $\nu$.

**Near incompressibility:** the rise in absolute stress between $\nu=0.40$ and $0.45$ may partly reflect
volumetric locking of P1 elements. No mixed formulation or $\nu$-specific convergence study was performed, so
the $\nu=0.45$ values are the least reliable.

---

## 7. Bellows geometry sensitivity

### 7.1 Opening angle (P; recomputed (F), strain column post-fix)

| $\Delta\theta$ | Max σ_e unreg. (Pa) | Max σ_e reg. (Pa) | Mean σ_e unreg. (Pa) | Mean σ_e reg. (Pa) | Max σ_e red. | Mean σ_e red. | Max strain red. |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.10° | 0.053962 | 0.021294 | 0.002353 | 0.001010 | 60.539 % | 57.079 % | 61.435 % |
| 0.25° | 0.134827 | 0.053201 | 0.005910 | 0.002582 | 60.541 % | 56.309 % | 61.438 % |
| 0.50° | 0.269396 | 0.106293 | 0.011936 | 0.005367 | 60.544 % | 55.037 % | 61.443 % |
| 1.00° | 0.537757 | 0.212147 | 0.024422 | 0.011584 | 60.550 % | 52.568 % | 61.453 % |

**Proportionality.** Going from $0.1^\circ$ to $1^\circ$, maximum stress grows by $9.97\times$ for both conditions,
nearly linear. Mean stress grows by $10.38\times$ (unregulated) and $11.47\times$ (regulated), faster than
linear. This is why the mean-stress reduction drifts from 57.1 % to 52.6 %. A plausible cause is examined in
§10.6.

### 7.2 Fold radius

Two variants of the study (METHODOLOGY §12.4): **original** (P), with thickness generated on $s/L_s$ (a smooth
ramp), and **corrected** (F), with thickness on $s$ in μm, consistent with all other experiments.

| $R_0$ (μm) | Max σ_e unreg. (orig.) | Max σ_e unreg. (corr.) | Max σ_e red. (orig.) | Max σ_e red. (corr.) | Mean σ_e red. (orig.) | Mean σ_e red. (corr.) |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 0.376198 | 0.463004 | 65.219 % | 66.742 % | 56.053 % | 55.878 % |
| 20 | 0.232342 | 0.269396 | 61.129 % | 60.544 % | 54.471 % | 55.037 % |
| 25 | 0.167634 | 0.191067 | 58.951 % | 57.389 % | 53.107 % | 54.041 % |
| 30 | 0.130917 | 0.147347 | 57.637 % | 55.193 % | 52.004 % | 53.190 % |
| 40 | 0.090988 | 0.101737 | 55.872 % | 52.892 % | 50.468 % | 51.854 % |

Corrected stress values were extracted from the final notebook's radius study **(D)**. The corrected run at
$R_0=20\ \mu$m reproduces the Section 10 baseline exactly (F).

**Both variants agree:** sharper folds produce higher stress, and thinning helps more on sharper folds. The
corrected variant shows a steeper decline of the max-stress benefit with radius (66.7 → 52.9 %, against
65.2 → 55.9 %). Conclusions about fold radius should use the corrected variant.

### 7.3 Thickness-to-radius scaling

Grouping the radius results by nominal $\chi=\bar h/R_0$, where pairs of cases share $\chi$:

| $\chi$ | Case A | Max σ_e (A) | Case B | Max σ_e (B) | Original A/B | Corrected A/B |
|---:|---|---:|---|---:|---:|---:|
| 0.333 | unreg., $R_0=30$ | orig. 0.130917 / corr. 0.147347 | reg., $R_0=15$ | orig. 0.130846 / corr. 0.153987 | 1.001 | 0.957 |
| 0.250 | unreg., $R_0=40$ | orig. 0.090988 / corr. 0.101737 | reg., $R_0=20$ | orig. 0.090313 / corr. 0.106293 | 1.007 | 0.957 |

With the smooth ramp profile, maximum stress depends almost only on $\chi$ (within 1 %). With realistic rough
profiles it collapses to within about 5 %, because local thickness at the fold, not the mean thickness, controls
the peak. $\chi$ is a useful organizing parameter, consistent with the bending-strain scaling
$\varepsilon\sim\chi\,\Delta\theta/(\pi-\theta_0)$, but it is not an exact similarity variable.

### 7.4 Cohesion-strength and strength-retention sensitivity (F)

Representative seed-42 pair, $60^\circ\rightarrow60.5^\circ$, $E=64$ Pa, $\nu=0.3$.

**At one assumed strength** ($\sigma_c=0.1$ Pa):

| | Max $FI$ | $f_{\mathrm{over}}$ | Overstressed area | Total area |
|---|---:|---:|---:|---:|
| Baseline thickness ($\alpha=1$) | 2.694 | 1.781 % | 76.2 μm² | 4278.3 μm² |
| Reduced thickness ($\alpha=0.5$) | 1.063 | 0.115 % | 2.5 μm² | 2174.3 μm² |

**Over assumed strengths 0.01–1 Pa** (200 values):

- $f_{\mathrm{over}}$ decreases as the assumed strength increases.
- The reduced-thickness biofilm is never more overstressed than the baseline: 0 of 200 strengths.
- Above $\sigma_c\approx0.27$ Pa, which is the baseline peak stress, neither biofilm exceeds the threshold.
- In particular, at Picioreanu et al.'s strengths of 0.4 and 1.0 Pa, there is no exceedance at this movement
  and stiffness. Because $FI\propto E\,\Delta\theta/\sigma_c$, this conclusion moves with the assumed $E$ and
  $\Delta\theta$.

**Strength retention.** For the reduced-thickness biofilm with strength $\beta\sigma_{c,\mathrm{baseline}}$, the
computed crossover $\beta^*$, below which thinning no longer lowers $f_{\mathrm{over}}$, is:

| Baseline $\sigma_c$ (Pa) | 0.03 | 0.05 | 0.10 | 0.20 | 0.25 |
|---|---:|---:|---:|---:|---:|
| Baseline $f_{\mathrm{over}}$ | 13.65 % | 9.11 % | 1.78 % | 0.26 % | 0.06 % |
| Crossover $\beta^*$ | 0.291 | 0.418 | 0.617 | 0.470 | 0.418 |

$\beta^*$ is not monotonic in the baseline strength and is not a universal constant. It depends on how the
stress is distributed over the element areas of the two biofilms. A visual reading of the colour map near
$\sigma_c=0.10$ Pa suggests a transition around $\beta\approx0.4$–$0.5$, because small positive $\Delta f$ is
nearly white. The computed crossover there is $\beta^*=0.617$, which the notebook marks on the map with the
$\Delta f=0$ contour. Where both fractions are zero (high strength, moderate $\beta$), $\Delta f=0$ means only
that neither biofilm exceeds the threshold. It does not mean the two are mechanically equivalent.

---

## 8. Interfacial traction convergence

Medium mesh (14 970 elements) vs. fine mesh (59 940 elements), seed 42 (P):

| Metric | Unreg. medium | Unreg. fine | δ unreg. | Reg. medium | Reg. fine | δ reg. |
|---|---:|---:|---:|---:|---:|---:|
| Max $t_n^+$ (Pa) | 0.044154 | 0.044789 | 1.42 % | 0.012854 | 0.013543 | 5.09 % |
| Max $\lvert t_s\rvert$ (Pa) | 0.041849 | 0.041775 | 0.18 % | 0.021182 | 0.020478 | 3.44 % |
| Mean $t_n^+$ | — | — | 1.99 % | — | — | 5.30 % |
| Mean $\lvert t_s\rvert$ | — | — | 0.02 % | — | — | 1.00 % |

| Reduction from thinning | Medium (F, matches P) | Fine (P) |
|---|---:|---:|
| Max $t_n^+$ | 70.89 % | 69.76 % |
| Max $\lvert t_s\rvert$ | 49.39 % | 50.98 % |
| Mean $t_n^+$ | 76.86 % | 76.05 % |
| Mean $\lvert t_s\rvert$ | 71.96 % | 71.69 % |

Production-mesh values (F): unregulated mean $t_n^+=0.005077$ Pa and mean $\lvert t_s\rvert=0.004841$ Pa;
regulated 0.001175 Pa and 0.001357 Pa.

**Reliability.**

- Interface tractions converge more slowly than bulk quantities, because they use the stress of a single
  boundary element, without extrapolation.
- Mean shear traction is the most reliable interface metric (≤ 1 %).
- Tensile normal traction on the thin regulated film is the least reliable (about 5 %).
- Peak values are more mesh-sensitive than length-weighted means.
- The *reductions* agree between the two meshes to within 1.6 percentage points, so the conclusion that thinning
  lowers interfacial loading is robust.
- Only two resolutions were compared, so no convergence order can be estimated.

---

## 9. FEM-to-grid conservation

$100\times100$ grid over the unregulated biofilm with 5 μm padding (F):

| Quantity | Relative error |
|---|---:|
| Total biofilm area $\sum\phi_{ij}\lvert G_{ij}\rvert$ vs. $\sum A_k$ | $0$ (round-off) |
| $\int\sigma_t\,dA$ | $2.4\times10^{-16}$ |
| $\int\sigma_e\,dA$ | $0$ (round-off) |

**What is not conserved or represented:**

- **Peaks and gradients.** Cell averaging smooths localized maxima, so the grid's maximum stress is lower than
  the FEM maximum.
- **Force and moment.** Integrating stress over area is not a force balance. Tractions are boundary quantities
  and are not mapped.
- **Tensor information.** Only scalar fields are mapped. Mapping tensor components requires a consistent
  rotation into the grid frame (see [INTEGRATION.md](INTEGRATION.md)).
- **Partial coverage.** If the grid does not cover the whole mesh, the area outside it is silently dropped, and
  the identities above no longer hold.

---

## 10. Scientific limitations

### 10.1 Biofilm viscoelasticity

Biofilms are viscoelastic. They store energy elastically on short timescales and flow or relax stress on long
ones.<sup>[L5–L8]</sup> Shaw et al. reported that elastic relaxation times are similar across a range of biofilms,
about 18 min.<sup>[L7]</sup> Klapper et al. described biofilm material as a viscoelastic fluid.<sup>[L8]</sup>

Our model is purely elastic. It represents the *instantaneous* response to a bellows movement that is fast
compared with relaxation. For slow movements, or for sustained deformation, stress would partly relax and the
elastic stresses here would overestimate it. The bellows movement rate in service has not been established for
this model.

### 10.2 Uncalibrated Young's modulus, $E=64$ Pa

This value comes from a 2D detachment model.<sup>[L1]</sup> That paper cites reported moduli of about
10–240 Pa. Later measurements span 70–700 Pa (OCT with fluid–structure modelling<sup>[L2]</sup>) and about
1–3 kPa (microfluidic rheometry<sup>[L3]</sup>); see also the review by Böl et al.<sup>[L4]</sup> Absolute
stresses scale with $E$ (§6.1), so they should be read as *relative* indicators, never compared with strength
values.

### 10.3 Uncalibrated Poisson's ratio, $\nu=0.3$

This value is also taken from [L1]. Hydrated gels are often closer to incompressible ($\nu\to0.5$), where the
P1 elements used here may lock (§6.2).

### 10.4 Homogeneous, isotropic material

Real biofilms are heterogeneous, with properties that vary with depth, density and matrix composition, and they
can be anisotropic.<sup>[L4]</sup> Local stiff or weak zones would redistribute stress.

### 10.5 Plane strain

The 2D section assumes the convolution is long and uniform out of plane, so there is no out-of-plane strain.
This does not hold near the ends of a convolution or for patchy biofilm coverage. Plane strain also produces a
nonzero $\sigma_{zz}=\nu(\sigma_{xx}+\sigma_{yy})$, which enters $\sigma_t$ and $\sigma_e$.

### 10.6 Small strain and finite rotation

Maximum principal strain reaches 0.90 % at $\Delta\theta=1^\circ$ and 0.45 % at the baseline $0.5^\circ$
(D, post-fix), within the small-strain range. Rotations are a separate issue. Each flank rotates rigidly by about $\Delta\theta/2$. Applied to a
finite rotation $\mathbf R(\omega)-\mathbf I$, linearized strain gives a spurious isotropic strain of about
$-\omega^2/2$. Under plane strain this produces a spurious von Mises stress of about $2\mu\,\omega^2/2$.

**Order-of-magnitude estimate (not a verified calculation):**

| $\Delta\theta$ | Spurious $\sigma_e$ | Compared with mean $\sigma_e$ |
|---|---:|---|
| 0.5° | ≈ $4.7\times10^{-4}$ Pa | about 4 % of unregulated, 9 % of regulated |
| 1° | ≈ $1.9\times10^{-3}$ Pa | about 8 % of unregulated, 16 % of regulated |

This is consistent with the faster-than-linear growth of *mean* stress, and with the drift of the mean-stress
reduction in §7.1. Peak stresses near the fold, which are dominated by genuine bending, are barely affected.
Larger bellows movements would require a geometrically nonlinear (finite-strain) formulation.

### 10.7 Idealized bellows geometry

A single symmetric convolution with straight flanks and a circular fold ($l=200$ μm, $R_0=20$ μm,
$\theta_0=60^\circ$) is an engineering idealization, not measured ISS hardware. The real corrugation profile,
pitch and wall thickness are unknown here. Curvature jumps at the tangent points are an artefact of the
idealization.

### 10.8 Prescribed substrate motion

The bellows deformation is imposed by an inextensibility rule (METHODOLOGY §4.4), not computed from metal
mechanics. Real convolutions may stretch, ovalize or deform asymmetrically. The biofilm does not influence the
metal; this is a one-way coupling, reasonable given the many orders of magnitude between metal and biofilm
stiffness.

### 10.9 Perfect adhesion

The interface can neither slip nor open. Interfacial tractions are computed, but there is no adhesion law, so
debonding cannot occur. With partial debonding, the internal stresses would be lower and the interface loads
redistributed.

### 10.10 No self-contact

Biofilm on the two flanks does not interact. Self-contact is irrelevant for the small opening changes used, but
closing the crevice substantially, or thicker films, could bring opposing surfaces into contact. The minimum
nonlocal surface separation was 4.6 μm in an earlier test geometry (P).

### 10.11 No fracture or detachment

No fracture mechanics, damage, cohesive-zone model or removal criterion is implemented. High $\sigma_e$,
$\sigma_t$ or $t_n^+$ marks **where loading is concentrated**, not where, when or how much biomass detaches.
The model makes no claim about fragment sizes or filter clogging.

### 10.12 Biological interventions may change material properties

The proposal's candidate thinning strategies (SinI/SinR control of matrix production, matrix-degrading enzymes,
γ-PGA-associated dispersal) act on the extracellular matrix, so they may change stiffness, cohesion and adhesion
as well as thickness. γ-PGA, for example, contributes to biofilm formation in some *B. subtilis*
strains.<sup>[L9, L10]</sup> We found no measured elastic modulus, cohesion or adhesion strength for biofilms
thinned by any of these routes. Assuming equal $E$ and $\nu$ for both members of a pair isolates the effect of
geometry. If an intervention softens the matrix, the
stresses in the regulated film would fall further, because stress is proportional to $E$; whether its
*strength* falls more than its stress is unknown. This "strength-retention" question was explored only as a
threshold-sensitivity map in the flat-substrate model (notebook 03), not for the bellows.

### 10.13 Hypothetical 50 % thickness reduction

$\alpha=0.5$ is a controlled geometric scenario, not a simulated or measured outcome of any genetic construct.
The FEM does not model gene expression or regulation. Uniform scaling also halves
roughness and slopes (METHODOLOGY §3.3). Other values of $\alpha$, or non-uniform thinning, were not studied.

### 10.14 Element-averaged stress recovery

P1 elements give piecewise-constant stress. Peaks and interface tractions are taken from single elements,
without smoothing or extrapolation, and depend on the mesh at the percent level (§4, §8). The reported maxima
include the strip ends, where the clamped base meets a free edge.

### 10.15 Limited experimental validation

No model output has been compared with measurements of our biofilm, of biofilm on bellows, or of ISS hardware
kinematics. Parameter values are illustrative and drawn from the literature on other organisms and systems.

---

### 10.16 Cohesive threshold analysis

- **No calibrated cohesion strength.** $\sigma_c$ and $\beta$ are scanned, not measured. No strength, nor its
  change under any proposed intervention, has been determined for our biofilm. Picioreanu et al.'s values
  (0.4, 1.0 Pa) were themselves assumptions for a different system.
- **Von Mises as a biofilm criterion.** The distortion-energy criterion comes from the yielding of ductile
  solids. It is insensitive to hydrostatic stress and does not distinguish tension from compression. Using it for
  a soft, viscoelastic, heterogeneous gel follows Picioreanu et al. but is unvalidated. Our 3D form also differs
  from their printed eq. (1); see METHODOLOGY §10.5.
- **No crack propagation or element deletion.** Overstressed elements are flagged, not removed. Stresses are
  not redistributed after "failure", so cascading failure, crack paths and post-failure geometry are not
  represented.
- **Cohesive versus adhesive failure.** The threshold is applied to bulk elements (cohesion). Adhesive failure
  at the metal is characterized only through interface tractions, for which no adhesion strength is known.
  Picioreanu et al. assumed the two strengths were equal; we make no such assumption, and we do not compare
  tractions with a strength.
- **Uncertain biological changes.** Matrix-targeting interventions may change stiffness as well as strength.
  The $\beta$ analysis varies strength only, at fixed $E$. Since $FI\propto E/\sigma_c$, a simultaneous change
  in $E$ would shift the crossover.
- **Convergence is not validity.** Mesh convergence of $\sigma_e$ (§4) supports the numerical stability of
  $f_{\mathrm{over}}$ near the stress peak. It does not show that exceeding the threshold corresponds to real
  failure.
- **Overstressed area is not detached biomass.** $f_{\mathrm{over}}$ is normalized by each biofilm's own area,
  so a larger fraction does not imply more vulnerable biomass in absolute terms. Neither the fraction nor the
  absolute overstressed area predicts detachment probability, detached mass or fragment size.

---

## 11. Defensible scientific conclusions

### What the model demonstrates

*These statements are true of the model as implemented, verified and converged.*

1. The FEM pipeline correctly imposes bellows-driven displacement on a curved biofilm. It is invariant to
   infinitesimal rigid motion, linear in $E$, and its grid mapping is conservative.
2. Under the stated assumptions, halving biofilm thickness reduces the maximum von Mises stress (mean 59.57 %),
   the area-weighted mean stress (52.22 %) and the maximum principal strain (60.37 %) in **all** 20 sampled
   morphologies.
3. In the representative case, maximum principal tensile stress falls by 62.22 %. Interfacial tensile normal
   traction falls by about 70–77 % and shear traction by about 49–72 % (peak and mean respectively). These
   reductions persist on a 4× finer mesh.
4. The relative reductions do not depend on $E$, shift by at most about 5 percentage points over $\nu=0.2$–$0.45$,
   and remain positive for every opening-angle change and fold radius tested.
5. At equal assumed cohesion strength (0.01–1 Pa), the reduced-thickness biofilm never has the larger
   overstressed-area fraction. With a strength-retention factor $\beta$ below a computed crossover $\beta^*$
   (0.29–0.62 for baseline strengths of 0.03–0.25 Pa), it does.

### What the model suggests

*These are plausible, but rest on unvalidated assumptions.*

1. Thinner biofilm attached to a flexing metal fold may experience substantially lower internal and interfacial
   mechanical loading, because a thinner layer is bent less for the same change in curvature (bending strain
   scales with $h/R_0$).
2. Sharper folds are mechanically more demanding locations for attached biofilm, and thickness reduction helps
   most there.
3. Reducing thickness lowers the stresses generated by bellows deformation, but whether this lowers cohesive
   vulnerability depends on the biofilm's strength. **Thickness and cohesion strength must be considered
   jointly.** Thinning helps if strength is largely retained, and its benefit can be offset by a sufficiently
   large loss of strength.

### What the model cannot establish

1. Whether any biofilm fractures, detaches or sloughs, how much biomass is released, or what size the fragments
   are. Threshold exceedance is not fragmentation.
2. Whether any proposed intervention (SinI/SinR matrix control, matrix-degrading enzymes, γ-PGA-associated
   dispersal) actually halves thickness, or how it alters stiffness, cohesion or adhesion.
3. Absolute stress levels in real ISS biofilms, given the uncalibrated $E$ and $\nu$ and the absence of
   viscoelasticity.
4. Behaviour under large bellows strokes, fluid loading, contact, or repeated cycling (fatigue).
5. Any effect on filter clogging or mission operations, or that thinning always reduces clogging risk.
6. The actual cohesion strength of our biofilm, or the value of $\beta$ produced by any intervention.

---

## References

Bibliographic details were checked against Crossref (title, authors, journal, volume, pages, DOI).

- **[L1]** Picioreanu C, van Loosdrecht MCM, Heijnen JJ. Two-dimensional model of biofilm detachment caused by internal stress from liquid flow. *Biotechnology and Bioengineering* 72(2):205–218 (2001). https://doi.org/10.1002/1097-0290(20000120)72:2<205::AID-BIT9>3.0.CO;2-L. Source of $E=64$ Pa and $\nu=0.3$; a local transcription is in `../docs/lit/sim/`.
- **[L2]** Picioreanu C, Blauert F, Horn H, Wagner M. Determination of mechanical properties of biofilms by modelling the deformation measured using optical coherence tomography. *Water Research* 145:588–598 (2018). https://doi.org/10.1016/j.watres.2018.08.070
- **[L3]** Hohne DN, Younger JG, Solomon MJ. Flexible microfluidic device for mechanical property characterization of soft viscoelastic solids such as bacterial biofilms. *Langmuir* 25(13):7743–7751 (2009). https://doi.org/10.1021/la803413x
- **[L4]** Böl M, Ehret AE, Bolea Albero A, Hellriegel J, Krull R. Recent advances in mechanical characterisation of biofilm and their significance for material modelling. *Critical Reviews in Biotechnology* 33(2):145–171 (2013). https://doi.org/10.3109/07388551.2012.679250
- **[L5]** Peterson BW, He Y, Ren Y, Zerdoum A, Libera MR, Sharma PK, et al. Viscoelasticity of biofilms and their recalcitrance to mechanical and chemical challenges. *FEMS Microbiology Reviews* 39(2):234–245 (2015). https://doi.org/10.1093/femsre/fuu008
- **[L6]** Gloag ES, Fabbri S, Wozniak DJ, Stoodley P. Biofilm mechanics: Implications in infection and survival. *Biofilm* 2:100017 (2020). https://doi.org/10.1016/j.bioflm.2019.100017
- **[L7]** Shaw T, Winston M, Rupp CJ, Klapper I, Stoodley P. Commonality of elastic relaxation times in biofilms. *Physical Review Letters* 93(9):098102 (2004). https://doi.org/10.1103/PhysRevLett.93.098102
- **[L8]** Klapper I, Rupp CJ, Cargo R, Purvedorj B, Stoodley P. Viscoelastic fluid description of bacterial biofilm material properties. *Biotechnology and Bioengineering* 80(3):289–296 (2002). https://doi.org/10.1002/bit.10376
- **[L9]** Stanley NR, Lazazzera BA. Defining the genetic differences between wild and domestic strains of *Bacillus subtilis* that affect poly-γ-DL-glutamic acid production and biofilm formation. *Molecular Microbiology* 57(4):1143–1158 (2005). https://doi.org/10.1111/j.1365-2958.2005.04746.x
- **[L10]** Yu Y, Yan F, Chen Y, Jin C, Guo J, Chai Y. Poly-γ-glutamic acids contribute to biofilm formation and plant root colonization in selected environmental isolates of *Bacillus subtilis*. *Frontiers in Microbiology* 7:1811 (2016). https://doi.org/10.3389/fmicb.2016.01811
- **[R1]** Gilbert et al. (2020), ICES-2020-566 (ISS wastewater tank operations; bellows-fold biomass accumulation). Cited from the team's modelling master plan; a local PDF is in `../docs/lit/mech/`. **Full bibliographic details not independently verified.**

No literature was found, or is cited, for the mechanical properties of *B. subtilis* biofilms thinned by any of the
proposed interventions,
for biofilm adhesion strength on stainless steel, or for ISS bellows kinematics.
