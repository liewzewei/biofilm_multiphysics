# Computational Methodology

**Structural mechanics of a biofilm on a deforming bellows convolution**

This document explains, from first principles, the model implemented in `src/` and demonstrated in
[`notebooks/FINAL_biofilm_mechanics.ipynb`](../notebooks/FINAL_biofilm_mechanics.ipynb). It describes what the
code *does*. How far the results can be trusted is covered in
[VALIDATION_AND_LIMITATIONS.md](VALIDATION_AND_LIMITATIONS.md), and coupling to the dispersion model in
[INTEGRATION.md](INTEGRATION.md).

Throughout, statements are tagged where it matters:

- **[Assumption]**: a choice about the physical system.
- **[Numerical]**: an approximation introduced by discretization or implementation.

---

## 1. Research objective

**Context.** Biofilm has accumulated in the wastewater tank of the International Space Station (ISS) Water
Processor Assembly. Wastewater fills the space between a metal bellows and the outer tank wall, and NASA
reports describe biomass collecting in the bellows folds.<sup>[R1]</sup> When the bellows expands or contracts,
its convolutions open and close, and biofilm stuck to the metal is deformed with it.

**Project context.** The iGEM project aims to engineer a non-pathogenic *Bacillus subtilis* biofilm that:

1. establishes a beneficial biofilm that competitively excludes undesirable microorganisms;
2. can be formed and thinned on demand through externally controlled regulation;
3. prevents excessive accumulation and so reduces the risk of downstream clogging.

These are the project's established objectives. The *means* of thinning is still a proposed implementation
strategy. The proposal's design iterations consider SinI/SinR control of extracellular-matrix production,
matrix-degrading enzymes and γ-PGA-associated dispersal, with external control by light and, at earlier stages,
by chemical induction. In some *B. subtilis* strains γ-PGA contributes to biofilm formation.<sup>[L9, L10]</sup>
None of these constructs has been experimentally validated in this project.

**Role of this model.** The FEM is the structural-mechanics component of the dry-lab modelling effort, alongside
models of flow and clogging, gene expression, ecology and competition, optogenetic kinetics, light penetrance
and γ-PGA-driven dispersion. It does **not** simulate gene expression or any regulatory mechanism. It simulates
the mechanical consequences of a *hypothetical change in biofilm thickness*: an *unregulated* (reference) profile
and a *regulated* (thinned) profile $h_r=\alpha h_u$ with $\alpha=0.5$. The thinned profile is a controlled
geometric assumption, not an established result of γ-PGA regulation or of any particular genetic circuit.
Biological interventions may alter mechanical properties as well as geometry (Section 3.3).

**Mechanical question.** Under identical bellows movement, how would a hypothetical halving of the biofilm
thickness change:

- the internal mechanical loading of the biofilm (von Mises stress, principal tensile stress and strain), and
- the loading transmitted across the biofilm–metal interface (normal and shear traction)?

We call these quantities *mechanical vulnerability indicators*. They measure **loading**, not fracture or
detachment, because no failure law is modelled.

---

## 2. Overview of the computational workflow

```
 random_height_profile          deform_rounded_bellows        build_curved_layer_mesh
 (stochastic thickness h(s)) ─┐ (60° → 60.5°, length kept) ─┐ (triangles on curved strip)
                              │                             │
                              └──► resample_bellows ────────┴──► solve_curved_biofilm ──► indicators
                                   (common material grid s)      (SfePy, plane strain)      │
                                                                                            ├─► calculate_interface_traction
                                                                                            └─► fem_to_grid (Cartesian grid)
```

| Stage | What it represents | Function (module) |
|---|---|---|
| Stochastic morphology | Biofilm thickness along the metal, $h(s)$ | `random_height_profile` (`src/stochastic_geometry.py`) |
| Bellows geometry | One rounded V-shaped convolution | `rounded_bellows_profile` (`src/rounded_bellows.py`) |
| Bellows kinematics | Opening by $\Delta\theta$ without stretching the metal | `deform_rounded_bellows`, `resample_bellows` (`src/rounded_bellows.py`) |
| Curved mesh | Biofilm strip offset along the surface normal | `build_curved_layer_mesh` (`src/curved_layer_mesh.py`) |
| Structural FEM | Plane-strain linear elasticity, prescribed base motion | `solve_curved_biofilm` (`src/curved_structural_fem.py`) |
| Indicators | σ_e, σ_t, ε_t, interface $t_n$ and $t_s$ | `solve_curved_biofilm`, `calculate_interface_traction` (`src/interface_traction.py`) |
| Grid integration | Element fields → Cartesian cells | `fem_to_grid` (`src/grid_mapping.py`) |

All geometry is specified in micrometres (μm). The FEM converts to SI units internally (Section 9).

> **Earlier flat-substrate model.** `src/layer_mesh.py` and `src/structural_fem.py` (`solve_biofilm`) implement
> the same elasticity problem for a biofilm $0\le y\le h(x)$ on a flat substrate, loaded by uniform surface
> traction or by prescribed horizontal substrate motion. That model was used in notebooks 02–04 and shares the
> constitutive and postprocessing equations of Section 7. Its results are summarized in
> [VALIDATION_AND_LIMITATIONS.md](VALIDATION_AND_LIMITATIONS.md). The rest of this document describes the
> curved bellows model.

---

## 3. Stochastic biofilm geometry

### 3.1 Material arc-length coordinate

Biofilm thickness is defined as a function of the **material arc-length coordinate** $s\in[0,L_s]$, measured
along the metal surface from one end of the convolution:

$$
h = h(s).
$$

**Why arc length rather than $x$?** A V-shaped surface is not a graph $y(x)$ over a single horizontal axis in
any useful sense, and its shape changes as the bellows opens. Arc length labels a *material point* on the
metal. Because the metal is assumed inextensible (Section 4.4), the same value of $s$ refers to the same metal
point in both the undeformed and the deformed configuration, so a thickness profile defined on $s$ stays
attached to the same material.

### 3.2 Random morphology generator

`random_height_profile(x, mean_height, roughness, n_features, min_height, rng)` builds a smooth random
fluctuation from $N$ Gaussian bumps:

$$
\tilde r(s)=\sum_{i=1}^{N} a_i\exp\!\left[-\frac12\left(\frac{s-\mu_i}{\sigma_i}\right)^2\right],
\qquad a_i\sim U(-1,1),\ \ \mu_i\sim U(0,L_s),\ \ \sigma_i\sim U(15,60).
$$

It then removes the mean and rescales to a prescribed RMS roughness $R_q$:

$$
r(s)=c\left[\tilde r(s)-\frac{1}{L_s}\int_0^{L_s}\tilde r\,ds\right],\qquad
c\ \text{chosen so that}\ \ \sqrt{\frac{1}{L_s}\int_0^{L_s} r^2\,ds}=R_q,
$$

$$
\boxed{h(s)=h_0+r(s)}.
$$

The integrals are evaluated with the trapezoidal rule on the sampling grid. A profile is accepted only if
$h(s)\ge h_{\min}$ everywhere; otherwise new random features are drawn (up to 1000 attempts).

**Why this construction?**

- **Gaussian features** give a smooth surface with a tunable lateral length scale, here 15–60 μm (an
  illustrative choice, not a measured morphology). The surface stays smooth enough to mesh without spurious
  corners.
- **Mean removal and RMS normalization** make every random profile carry *exactly* the same biomass
  ($h_0L_s$) and the same roughness amplitude. Differences between samples are then due only to the
  arrangement of features, which is what the ensemble in Section 11 is designed to probe.
- **The minimum-thickness constraint** prevents zero-thickness gaps, which would disconnect the mesh.

**Important unit convention.** The feature widths $\sigma_i\in[15,60]$ are drawn in the **same units as the
input coordinate**. With $s$ in μm they are physical widths in micrometres. See Section 12.4 for the
consequence of passing a normalized coordinate.

**Parameters used** (bellows model): $h_0=10\ \mu$m, $R_q=3\ \mu$m, $N=10$, $h_{\min}=2\ \mu$m.
The representative morphology uses seed 42; the ensemble uses seeds 0–19 via `np.random.default_rng(seed)`.

### 3.3 Paired thickness scenarios (unregulated / regulated)

For every random reference (unregulated) profile $h_u(s)$ we construct a thinned (regulated) counterpart

$$
\boxed{h_r(s)=\alpha\,h_u(s),\qquad \alpha=0.5}.
$$

**[Assumption] Uniform thickness scaling.** The thinned scenario shrinks the biofilm proportionally
everywhere. This is a controlled geometric scenario, not a model of how any biological intervention actually
thins a biofilm. It keeps the *shape* of the morphology, so each pair differs only by the thickness change, a
controlled comparison. As a consequence, absolute roughness ($R_q\to\alpha R_q$) and surface slopes also scale
by $\alpha$. The experiment therefore tests "the same biofilm, uniformly thinner", not "a different biofilm with
lower mean thickness".

**[Assumption] Unchanged material.** Both members of a pair have the same $E$ and $\nu$. A real thinning
intervention, such as SinR-mediated matrix repression, matrix-degrading enzymes or γ-PGA-associated dispersal,
could also change matrix stiffness, cohesion and adhesion; that possibility is outside the present model (see
[VALIDATION_AND_LIMITATIONS.md §10](VALIDATION_AND_LIMITATIONS.md#10-scientific-limitations)).

---

## 4. Bellows geometry and kinematics

### 4.1 Rounded V-shaped convolution

**[Assumption]** One convolution is idealized as two straight segments of length $l$ meeting at opening angle
$\theta$, joined by a circular arc (the *fold*) of radius $R$ that is tangent to both segments.

`rounded_bellows_profile(segment_length, angle_deg, radius, n_straight=100, n_arc=101)` places the fold's centre
at the origin. Each straight segment makes the angle

$$
\varphi=\frac{\pi-\theta}{2}
$$

with the horizontal, and the arc is parameterized by

$$
\mathbf r_{\mathrm{arc}}(\psi)=R\begin{bmatrix}\cos\psi\\ \sin\psi\end{bmatrix},\qquad
\psi\in\left[\tfrac{\pi}{2}-\varphi,\ \tfrac{\pi}{2}+\varphi\right],
$$

so that its endpoints, the tangent points, join the straight segments with a continuous tangent. The profile is
sampled with 100 points per segment and 101 points on the arc (300 points in total). It is translated so that
its lowest point is at $y=0$ and returned together with the cumulative chord length $s$, a discrete
approximation of arc length.

### 4.2 Why a rounded fold?

A sharp V would have infinite curvature at the apex. Any biofilm attached there would see a geometric stress
singularity whose peak value depends on the mesh rather than on the physics. A finite radius $R$ makes the
curvature finite and turns the fold radius into a meaningful design parameter (Section 12.4).

### 4.3 Curvature

Curvature measures how quickly the tangent direction turns per unit arc length:

$$
\kappa=\left\lVert\frac{d\mathbf t}{ds}\right\rVert,\qquad
\kappa(s)=\begin{cases}0 & \text{on the straight segments}\\ 1/R & \text{on the fold}\end{cases}
$$

The tangent is continuous everywhere, but $\kappa$ jumps from $0$ to $1/R$ at the two tangent points. Opening
the bellows changes $\kappa$ only on the fold. That is where an attached layer must *bend*, while biofilm on
the flanks is mostly carried along rigidly.

### 4.4 Material-length preservation and the radius law

**[Assumption]** The metal is inextensible: as the bellows opens from $\theta_0$ to $\theta_1$, the straight
segments keep length $l$ and the fold keeps its arc length.

Between the two flanks the tangent direction turns through the angle $\pi-\theta$ (the supplement of the
opening angle), so the fold arc length is

$$
L_{\mathrm{fold}}=R\,(\pi-\theta).
$$

Holding $L_{\mathrm{fold}}$ fixed between the two configurations gives

$$
R_0(\pi-\theta_0)=R_1(\pi-\theta_1)
\quad\Longrightarrow\quad
\boxed{R_1=R_0\,\frac{\pi-\theta_0}{\pi-\theta_1}}
$$

This relation is implemented in `deform_rounded_bellows`. The total length
$L_s=2l+R(\pi-\theta)$ is then the same in both configurations. Opening the bellows ($\theta_1>\theta_0$)
reduces the turning angle, so the fold radius increases and the fold flattens. The exact change in fold
curvature is

$$
\Delta\kappa=\frac{1}{R_1}-\frac{1}{R_0}=\frac{(\pi-\theta_1)-(\pi-\theta_0)}{R_0(\pi-\theta_0)}
=-\frac{\theta_1-\theta_0}{R_0(\pi-\theta_0)}.
$$

For the baseline $R_0=20\ \mu$m and $60^\circ\rightarrow60.5^\circ$: $R_1=20.0837\ \mu$m and
$\Delta\kappa\approx-2.08\times10^{-4}\ \mu\mathrm m^{-1}$.

**Why prescribe the kinematics geometrically?** We do not know the bellows' wall thickness, alloy or
corrugation profile, so computing its deformation from metal mechanics is out of reach. Length preservation is
the simplest kinematic rule consistent with a thin metal sheet that bends far more easily than it stretches.

### 4.5 Material-point correspondence and the substrate displacement

Because both configurations have the same total length, the material point labelled $s$ in the undeformed
bellows is the point at arc length $s$ in the deformed bellows. `resample_bellows(points, s, s_new)` linearly
interpolates each configuration onto the *same* uniform grid $s_j=j\,L_s/(N-1)$, with $N=500$. The substrate
displacement of every material point is then

$$
\boxed{\mathbf u_b(s_j)=\mathbf r_b(s_j;\theta_1)-\mathbf r_b(s_j;\theta_0)}.
$$

For the baseline movement, $\max\lvert\mathbf u_b\rvert=0.783\ \mu$m.

**Reference frame.** Both profiles are generated symmetric about $x=0$ with their free ends at the same
height. This fixes the rigid-body part of $\mathbf u_b$. Any other choice differs by a rigid motion, which
produces no strain in linear elasticity to first order (verified in
[VALIDATION_AND_LIMITATIONS.md §3](VALIDATION_AND_LIMITATIONS.md#3-fem-verification)).

### 4.6 The upward-opening crevice

The generator produces an inverted V (fold at the top). Biofilm in the bellows grows on the inner, concave side
of the convolution, which we represent as an upward-opening V by a vertical reflection:

$$
x_{\mathrm{new}}=x,\qquad y_{\mathrm{new}}=y^{(0)}_{\max}-y.
$$

Both configurations are reflected about the **same** level $y^{(0)}_{\max}$ and shifted by the **same**
constant. This preserves distances, arc lengths and $\mathbf u_b(s)$, because the reflection is an isometry and
the shift is common to both configurations. In the notebook this is the helper `crevice_configurations`, which
only composes the functions above.

---

## 5. Curved biofilm geometry

### 5.1 Tangent and normal

For the substrate curve $\mathbf r_b(s)$, the unit tangent and the left-hand unit normal are

$$
\mathbf t(s)=\frac{d\mathbf r_b/ds}{\lVert d\mathbf r_b/ds\rVert},\qquad
\mathbf n(s)=\begin{bmatrix}-t_y(s)\\ t_x(s)\end{bmatrix}.
$$

`build_curved_layer_mesh` estimates $d\mathbf r_b$ with central differences (`np.gradient`) on the sampled
points and normalizes the result. With the substrate ordered from left to right, the left-hand normal
(`normal_side="left"`) points *into* the crevice, where the biofilm grows. `normal_side="right"` flips it.

### 5.2 Curvilinear coordinates and the normal-offset construction

Every point of the biofilm is located by its arc-length position $s$ and a normalized through-thickness
coordinate $\eta\in[0,1]$:

$$
\boxed{\mathbf r(s,\eta)=\mathbf r_b(s)+\eta\,h(s)\,\mathbf n(s)}.
$$

- $\eta=0$ is the biofilm–metal interface $\Gamma_s$.
- $\eta=1$ is the exposed biofilm surface $\Gamma_{\mathrm{top}}$.
- $s=0$ and $s=L_s$ are the two lateral ends of the strip.

**Why this construction?** It is the natural generalization of the flat-substrate layer $0\le y\le h(x)$:
thickness is measured perpendicular to the metal, which is how a coating grows. It also yields a structured
$(s,\eta)$ grid that maps directly to a mesh (Section 6), and keeps the substrate nodes exactly on the
bellows points where $\mathbf u_b$ is known.

**Validity condition.** On a concave surface, offsetting by more than the local radius of curvature folds the
offset curve over itself. The construction is therefore valid only if

$$
h(s)\,\kappa(s)<1
$$

locally, and if the exposed surfaces of the two flanks do not meet. In the baseline, $\max h=15.5\ \mu$m is
less than $R_0=20\ \mu$m. These conditions are not enforced by the mesher; they are checked after meshing
(Section 6.4).

---

## 6. Mesh generation

### 6.1 Boundary-fitted triangular mesh

`build_curved_layer_mesh(substrate_points, thickness, ny, normal_side)` samples $\eta$ uniformly,
$\eta_k=k/n_y$, at each of the $N$ substrate points. This creates $N(n_y+1)$ nodes. The node for
$(s_i,\eta_k)$ has index $k\,N+i$, so layer $k=0$ (indices $0\ldots N-1$) is exactly the substrate.

### 6.2 Element connectivity and orientation

Each quadrilateral cell $(i,k)$ with corners $n_{00},n_{10},n_{01},n_{11}$ is split along the same diagonal
into the triangles $[n_{00},n_{10},n_{11}]$ and $[n_{00},n_{11},n_{01}]$. This gives $2(N-1)n_y$ linear
triangles. Each triangle's signed area

$$
A_e=\tfrac12\left[(x_2-x_1)(y_3-y_1)-(y_2-y_1)(x_3-x_1)\right]
$$

is computed, and any clockwise triangle ($A_e<0$) has two vertices swapped. **Why:** FEM assembly assumes a
consistent counter-clockwise orientation, so that element Jacobians are positive.

The function also returns boundary node sets: `"substrate"`, `"top"`, `"left"` and `"right"`.

Production mesh: $N=500$ and $n_y=15$, giving 8 000 nodes and 14 970 triangles per biofilm.

### 6.3 Element-quality metrics

`evaluate_mesh_quality(points, triangles)` (`src/mesh_quality.py`) returns:

- the signed area of each element (positivity means no inversion);
- the minimum interior angle, computed with the law of cosines. Small angles degrade conditioning and
  stress-recovery accuracy;
- the aspect ratio, longest edge divided by shortest edge.

### 6.4 Uniform arc-length resampling

The raw bellows profile is sampled more densely on the arc than on the straight segments (spacing 0.42–2.02 μm).
Meshing it directly produced minimum element angles of about 7°. Resampling onto a uniform $s$-grid
(spacing 0.886 μm) raised the minimum angle to 20.6° for the unregulated mesh ($n_y=15$). The regulated mesh,
being half as thick, has a minimum angle of 10.7°. **[Numerical]** Resampling replaces the arc with a
500-point polyline. The total length changes by less than $10^{-7}$ relative.

### 6.5 Geometric validity checks and their limits

After meshing, the notebook checks that:

1. all signed areas are positive;
2. the boundary polygon (substrate followed by the reversed exposed surface) is valid according to Shapely;
3. the exposed surface is a simple curve, with no self-intersection;
4. the summed element area equals the polygon area (to about $10^{-15}$ relative).

These checks catch inverted elements and a folded exposed surface. They do **not** prove that no two interior
elements overlap, and the mesher does not enforce $h\kappa<1$ itself.

---

## 7. Governing mechanical equations

### 7.1 Quasi-static equilibrium

**[Assumption]** Bellows motion is slow compared with elastic wave propagation in the biofilm, and body forces
are negligible: there is no gravity in orbit, and fluid pressure is not modelled. Linear momentum balance then
reduces to

$$
\nabla\cdot\boldsymbol\sigma=\mathbf 0\quad\text{in }\Omega.
$$

### 7.2 Small-strain kinematics

**[Assumption]** Displacement gradients are small, $\lVert\nabla\mathbf u\rVert\ll1$, so strain is the
symmetric part of the displacement gradient:

$$
\boldsymbol\varepsilon=\tfrac12\left(\nabla\mathbf u+\nabla\mathbf u^{T}\right)
=\begin{bmatrix}\varepsilon_{xx}&\varepsilon_{xy}\\ \varepsilon_{xy}&\varepsilon_{yy}\end{bmatrix},
\qquad \gamma_{xy}=2\varepsilon_{xy}.
$$

**Why:** this makes the problem linear. The antisymmetric part of $\nabla\mathbf u$ (infinitesimal rotation)
produces no strain. That is exact only for *infinitesimal* rotations, which is why bellows movements are kept
small (Section 12.3).

### 7.3 Linear isotropic elasticity

**[Assumption]** The biofilm is homogeneous, isotropic and linearly elastic. Stress depends linearly on
strain through two constants, Young's modulus $E$ (stiffness) and Poisson's ratio $\nu$ (lateral contraction):

$$
\boldsymbol\sigma=\lambda\,\mathrm{tr}(\boldsymbol\varepsilon)\,\mathbf I+2\mu\,\boldsymbol\varepsilon,
\qquad
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\quad \mu=\frac{E}{2(1+\nu)}.
$$

### 7.4 Plane strain

**[Assumption]** The modelled cross-section is part of a structure that is long in the out-of-plane direction
$z$ (along a bellows convolution), so $\varepsilon_{zz}=\varepsilon_{xz}=\varepsilon_{yz}=0$. In Voigt
notation with engineering shear strain:

$$
\begin{bmatrix}\sigma_{xx}\\ \sigma_{yy}\\ \tau_{xy}\end{bmatrix}
=\underbrace{\frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}1-\nu&\nu&0\\ \nu&1-\nu&0\\ 0&0&\tfrac{1-2\nu}{2}\end{bmatrix}}_{\mathbf D}
\begin{bmatrix}\varepsilon_{xx}\\ \varepsilon_{yy}\\ \gamma_{xy}\end{bmatrix},
\qquad
\boxed{\sigma_{zz}=\nu(\sigma_{xx}+\sigma_{yy})}.
$$

Plane strain removes out-of-plane **strain**, not out-of-plane **stress**. $\sigma_{zz}$ is generally nonzero and
is included in the von Mises and principal-stress indicators. The matrix $\mathbf D$ is built by
`stiffness_from_youngpoisson(2, E, nu, plane="strain")`. For $E=64$ Pa and $\nu=0.3$:

$$
\mathbf D=\begin{bmatrix}86.15&36.92&0\\36.92&86.15&0\\0&0&24.62\end{bmatrix}\ \mathrm{Pa}.
$$

### 7.5 Material parameters

| Parameter | Value | Status |
|---|---|---|
| $E$ | 64 Pa | Illustrative; the value used by Picioreanu et al.<sup>[L1]</sup> in a 2D biofilm detachment model |
| $\nu$ | 0.3 | Illustrative; same source |

Measured biofilm moduli span orders of magnitude, from tens of Pa to kPa.<sup>[L1–L4]</sup> Under prescribed
displacement, stresses scale exactly with $E$, while strains and percentage reductions do not depend on it
(Section 12.1). The choice of $E$ therefore affects absolute stresses only.

---

## 8. Boundary conditions

| Boundary | Condition | Physical meaning |
|---|---|---|
| Interface $\Gamma_s$ ($\eta=0$) | $\mathbf u(s,0)=\mathbf u_b(s)$, both components | **[Assumption]** Perfect adhesion: no slip or separation at the metal |
| Exposed surface $\Gamma_{\mathrm{top}}$ ($\eta=1$) | $\boldsymbol\sigma\mathbf n=\mathbf 0$ | No fluid traction; bellows motion is the only load |
| Lateral ends ($s=0,\ L_s$) | $\boldsymbol\sigma\mathbf n=\mathbf 0$ | Free ends of the modelled strip |

**Reference configuration.** The biofilm is stress-free in the undeformed crevice ($\theta_0=60^\circ$), and the
mesh is built in that configuration. The FEM solves for the displacement *from* this state when the base moves
by $\mathbf u_b$.

**Material-point correspondence.** The substrate nodes of the mesh *are* the resampled reference bellows points
$\mathbf r_b(s_j;\theta_0)$. Each node therefore receives the displacement of its own material point. In
`solve_curved_biofilm`, the class `PrescribedBellowsDisplacement` assigns $\mathbf u_b$ to each boundary node by
nearest-neighbour lookup. It raises an error if any node lies more than $10^{-10}$ m from a substrate point, so
a mismatched geometry cannot be silently interpolated.

**Why Dirichlet rather than traction loading?** The bellows motion is a kinematic input, and the metal is far
stiffer than the biofilm. Prescribing displacement is equivalent to assuming the biofilm cannot push back on
the metal.

---

## 9. FEM implementation

### 9.1 Framework and units

`solve_curved_biofilm(points, triangles, boundaries, substrate_displacement, E=64.0, nu=0.3)` uses
**SfePy**. Inputs are in μm. The function multiplies coordinates and prescribed displacements by `UM = 1e-6`
before building the SfePy mesh, so assembly, solution and outputs are in SI units: m for displacement, Pa
for stress, m² for element areas.

### 9.2 Weak formulation

Multiplying equilibrium by a test displacement $\mathbf v$ that vanishes on $\Gamma_s$, integrating by parts, and
using the traction-free conditions on all other boundaries gives the principle of virtual work. Find $\mathbf u$
with $\mathbf u=\mathbf u_b$ on $\Gamma_s$ such that

$$
\int_\Omega\boldsymbol\varepsilon(\mathbf v)^{T}\,\mathbf D\,\boldsymbol\varepsilon(\mathbf u)\,d\Omega=0
\qquad\forall\,\mathbf v,\ \ \mathbf v|_{\Gamma_s}=\mathbf 0 .
$$

The right-hand side is zero because there are no applied tractions or body forces. The only driving term is
the prescribed boundary displacement. In SfePy this is the term `dw_lin_elastic(biofilm.D, v, u)` on the
region `Omega`.

### 9.3 Linear triangular elements

**[Numerical]** The displacement field is approximated with first-order Lagrange (P1) shape functions on each
triangle:

$$
\mathbf u(\mathbf x)\approx\sum_{a=1}^{3}N_a(\mathbf x)\,\mathbf u_a ,
$$

so displacement is continuous and piecewise linear, and **strain and stress are constant within each
element**. The element stiffness $\mathbf K_e=\int_{\Omega_e}\mathbf B^{T}\mathbf D\,\mathbf B\,d\Omega$ is
integrated with a second-order quadrature rule (`Integral("i", order=2)`), which is exact for P1.

### 9.4 Assembly and solution

Global assembly gives $\mathbf K\mathbf U=\mathbf F$. After eliminating the 1 000 Dirichlet degrees of
freedom (500 substrate nodes × 2), the production mesh has 15 000 unknowns. The system is linear and is solved
by a single Newton step (`Newton(i_max=1, eps_a=1e-10)`) with a sparse direct solver (`ScipyDirect`). The
nonlinear residual falls by about 15 orders of magnitude in that step (relative residual
$\approx2.5\times10^{-15}$), confirming an exact linear solve.

### 9.5 Displacement extraction

The nodal displacement vector is reshaped to `(N_nodes, 2)` in metres (`result["displacement"]`), and its
magnitude $\lvert\mathbf u\rvert$ is stored as `result["displacement_magnitude"]`. The maximum, in μm, is
`result["max_displacement_um"]`.

### 9.6 Stress and strain postprocessing

- **Strain:** `ev_cauchy_strain` in `el_avg` mode gives one strain vector per element, stored as
  `result["strain"]` with shape `(M, 3)`. **SfePy stores the shear component doubled**, so column 2 is the
  engineering strain $\gamma_{xy}=2\varepsilon_{xy}$. The solvers halve it before computing principal strains
  (Section 10.3).
- **Stress:** `ev_cauchy_stress` in `el_avg` mode gives $[\sigma_{xx},\sigma_{yy},\tau_{xy}]$ per element, stored
  as `result["stress"]` with shape `(M, 3)` in Pa. $\sigma_{zz}=\nu(\sigma_{xx}+\sigma_{yy})$ is computed from these.
- **Areas:** element areas in m² (`result["areas"]`), used for area-weighted means.

**[Numerical]** All derived indicators are element averages. P1 stresses are discontinuous across elements, and
no smoothing or superconvergent recovery is applied.

---

## 10. Mechanical indicators

All indicators are computed per element inside `solve_curved_biofilm`, except the interface tractions. The
symbols here match the code.

### 10.1 Von Mises equivalent stress

$$
\sigma_e=\sqrt{\tfrac12\left[(\sigma_{xx}-\sigma_{yy})^2+(\sigma_{yy}-\sigma_{zz})^2+(\sigma_{zz}-\sigma_{xx})^2\right]+3\tau_{xy}^2}
$$

Output: `sigma_e` per element, and `global_peak`, its maximum over all elements.

**Meaning:** a single scalar measuring the shear (distortional) part of the stress state. It ignores hydrostatic
pressure, and its sign is lost.

**Limitation:** the von Mises criterion was developed for yielding of ductile metals. For a biofilm it is used
only as a comparative loading indicator, never as a failure criterion.

### 10.2 Principal tensile stress

The in-plane principal stresses are the eigenvalues of the 2D stress tensor:

$$
\sigma_{1,2}=\frac{\sigma_{xx}+\sigma_{yy}}{2}\pm\sqrt{\left(\frac{\sigma_{xx}-\sigma_{yy}}{2}\right)^2+\tau_{xy}^2}.
$$

Under plane strain, $\sigma_{zz}$ is the third principal stress, so the maximum principal tensile stress is

$$
\boxed{\sigma_t=\max\left(0,\ \sigma_1,\ \sigma_2,\ \sigma_{zz}\right)}
$$

Output: `principal_tensile_stress` per element, and `max_tensile_stress`.

**Meaning:** soft cohesive materials tend to tear under tension rather than compression. $\sigma_t$ keeps only
tensile loading and discards compression, which von Mises would count equally.

**Limitation:** no tensile strength of our biofilm is known, so $\sigma_t$ cannot be converted into a failure
probability.

### 10.3 Principal tensile strain

The principal strains are the eigenvalues of the strain tensor:

$$
\varepsilon_{1,2}=\frac{\varepsilon_{xx}+\varepsilon_{yy}}{2}\pm\sqrt{\left(\frac{\varepsilon_{xx}-\varepsilon_{yy}}{2}\right)^2+\varepsilon_{xy}^2},
\qquad
\varepsilon_t=\max(0,\varepsilon_1),
$$

with $\varepsilon_{zz}=0$ under plane strain. The code also reports
$\max_e\max(\lvert\varepsilon_1\rvert,\lvert\varepsilon_2\rvert)$ as `max_principal_strain`, which serves as
the small-strain diagnostic. Outputs: `principal_tensile_strain`, `max_tensile_strain` and
`max_principal_strain`.

> **Implementation note.** Because SfePy returns $\gamma_{xy}$, both `solve_curved_biofilm` and
> `solve_biofilm` use $\varepsilon_{xy}=\gamma_{xy}/2$ in the formula above. Before Stage 4 the solvers used
> $\gamma_{xy}$ directly, which overestimated principal strains. That error was fixed, and
> `tests/test_strain_convention.py` now checks that principal strains are invariant under rotation (see
> [VALIDATION_AND_LIMITATIONS.md §3.5](VALIDATION_AND_LIMITATIONS.md#35-known-inconsistency-shear-strain-convention)).

### 10.4 Interfacial normal and shear traction

`calculate_interface_traction(points, triangles, stress, substrate_nodes)` evaluates Cauchy's formula on each
substrate edge $e$ between consecutive substrate nodes:

$$
\mathbf t=\boldsymbol\sigma\,\mathbf n,\qquad
\boxed{t_n=\mathbf n^{T}\boldsymbol\sigma\,\mathbf n},\qquad
\boxed{t_s=\mathbf s^{T}\boldsymbol\sigma\,\mathbf n},
$$

where $\mathbf s$ is the unit edge tangent and $\mathbf n=(-s_y,s_x)$ its left-hand normal. $\boldsymbol\sigma$ is
the in-plane stress tensor of the **single** element adjacent to that edge; the function raises an error if the
edge has zero or two neighbouring elements. Results are returned per edge: `midpoints`, `lengths`,
`normal_traction` and `shear_traction`.

The reported adhesion-loading indicators are:

$$
t_n^{+}=\max(t_n,0)\ \ \text{(tensile, pulling the biofilm off the metal)},\qquad \lvert t_s\rvert\ \ \text{(interfacial shear)},
$$

together with their maxima and length-weighted means
$\overline{q}=\sum_e\ell_e q_e/\sum_e\ell_e$.

**Orientation.** The left-hand normal points into the biofilm only for meshes built with
`normal_side="left"`, which every simulation uses. The final notebook verifies this explicitly. $t_n$ is
unchanged if $\mathbf n$ is reversed, and $\lvert t_s\rvert$ is unchanged if $\mathbf s$ is reversed.

**Limitations:** element-averaged stresses are evaluated at the interface without extrapolation, so peak
tractions depend on the mesh (see [VALIDATION_AND_LIMITATIONS.md §8](VALIDATION_AND_LIMITATIONS.md#8-interfacial-traction-convergence)).
No adhesion strength is known, so these values indicate where and how strongly the interface is loaded, not
whether it debonds.

---

## 11. Stochastic comparison

### 11.1 Paired design

For each seed $k=0,\dots,19$:

1. generate $h_u^{(k)}(s)$ (Section 3.2) and set $h_r^{(k)}=0.5\,h_u^{(k)}$;
2. mesh both on the *same* reference crevice;
3. apply the *same* $\mathbf u_b(s)$ for $60^\circ\rightarrow60.5^\circ$, with $E=64$ Pa and $\nu=0.3$;
4. record the maximum principal strain, the maximum von Mises stress and the area-weighted mean von Mises stress.

**Why paired?** Morphology alone changes the stresses, because peaks concentrate under thick features near the
fold. Comparing each regulated biofilm with *its own* unregulated parent removes this morphological variance
from the comparison, so the spread of paired reductions measures how consistent the thinning effect is.

### 11.2 Controlled loading

Loading is identical within each pair because the base displacement is a kinematic input independent of
biofilm thickness. No traction normalization is needed, unlike the earlier flat-substrate traction study, which
normalized by surface length.

### 11.3 Percentage reduction metric

$$
\mathrm{reduction}_Q=100\left(1-\frac{Q_r}{Q_u}\right)\% .
$$

Notebook 04 reported the equivalent negative "change" $100(Q_r/Q_u-1)\%$. A positive reduction means the
regulated biofilm has the smaller value. Ensemble statistics (mean, median, min, max, number of pairs reduced)
are computed across the 20 per-pair values.

### 11.4 Interpretation

A consistent sign and narrow spread across morphologies means the effect is a property of the thickness change
under these assumptions, not of one particular random shape. The ensemble does **not** sample uncertainty in
material parameters, bellows geometry or the degree of thinning $\alpha$; the first two are explored one at a
time in Section 12, and $\alpha$ is fixed at 0.5.
The results are reported in [VALIDATION_AND_LIMITATIONS.md §5](VALIDATION_AND_LIMITATIONS.md#5-stochastic-robustness).

---

## 12. Sensitivity analysis

Each study varies one parameter about the baseline ($R_0=20\ \mu$m, $60^\circ\rightarrow60.5^\circ$, $E=64$ Pa,
$\nu=0.3$, seed 42) and compares regulated with unregulated biofilms.

### 12.1 Young's modulus

With displacement prescribed on $\Gamma_s$ and no other loads, the weak form (Section 9.2) is homogeneous in
$\mathbf D\propto E$. Scaling $E$ therefore leaves $\mathbf u$ and $\boldsymbol\varepsilon$ unchanged and scales
$\boldsymbol\sigma$ exactly by the same factor:

$$
\mathbf u(E)=\mathbf u(E_0),\qquad \boldsymbol\sigma(E)=\frac{E}{E_0}\,\boldsymbol\sigma(E_0).
$$

Percentage reductions are therefore independent of $E$, as long as regulated and unregulated biofilms share the
same $E$. The study ($E=32$–$6400$ Pa) serves as a check of this linearity.

### 12.2 Poisson's ratio

$\nu$ controls the coupling between normal strain components and the magnitude of $\sigma_{zz}$, so it changes
the stress distribution and not just its scale. Values $\nu\in\{0.2,0.3,0.4,0.45\}$ were tested.

**[Numerical]** As $\nu\to0.5$ the material becomes incompressible and P1 displacement elements can *lock*
(become artificially stiff). Results at $\nu=0.45$ should be read with caution.

### 12.3 Bellows opening angle

$\Delta\theta\in\{0.1^\circ,0.25^\circ,0.5^\circ,1^\circ\}$ starting from $60^\circ$. For a linear model the stress
should scale with the size of the prescribed displacement. Departures from proportionality indicate that
$\mathbf u_b(\Delta\theta)$ itself is not linear in $\Delta\theta$, for example because of finite flank rotation.
This is why the baseline uses a small movement of $0.5^\circ$.

### 12.4 Fold radius

$R_0\in\{15,20,25,30,40\}\ \mu$m. For each radius, a new bellows geometry, arc-length grid ($N=500$) and
displacement field are generated with the same seed-42 thickness parameters.

Two variants exist:

- **Original (notebook 04):** thickness was generated on the *normalized* coordinate $s/L_s\in[0,1]$. Because the
  feature widths are drawn in input units (15–60), every Gaussian is far wider than the domain, and the profile
  becomes a smooth monotonic ramp (about 4→15 μm) rather than a rough morphology.
- **Corrected (final notebook, Section 15b):** thickness was generated on $s$ in μm, consistent with every other
  experiment. At $R_0=20\ \mu$m this reproduces the baseline exactly. Because $L_s$ changes with $R_0$, the
  corrected profiles differ slightly from one radius to the next.

### 12.5 Thickness-to-radius ratio

Thin-layer bending theory suggests that a layer of thickness $h$ on a surface whose curvature changes by
$\Delta\kappa$ experiences bending strains of order $z\,\Delta\kappa$, where $z\le h$ is the distance from the
interface. Using the exact $\Delta\kappa$ from Section 4.4:

$$
\lvert\varepsilon_{\mathrm{bend}}\rvert\sim h\,\lvert\Delta\kappa\rvert=\frac{h}{R_0}\,\frac{\lvert\Delta\theta\rvert}{\pi-\theta_0}.
$$

This motivates plotting the stress against $\chi=\bar h/R_0$ (nominal $\bar h=10$ μm or $5$ μm). It is a
**scaling argument only**. It ignores the clamped interface, the curvature jumps at the tangent points and the
2D stress state. For the baseline it predicts $\varepsilon\sim2\times10^{-3}$, the same order of magnitude as
the FEM peak strains.

---

## 13. FEM-to-grid integration

**Purpose.** The dispersion model runs on a Cartesian grid. `fem_to_grid` transfers element-averaged FEM fields
onto that grid without changing either model. Full API details are in [INTEGRATION.md](INTEGRATION.md).

### 13.1 Cartesian grid projection

For each grid cell $G_{ij}$ (row $j$ along $y$, column $i$ along $x$) and each FEM triangle $T_k$, the exact
intersection area $\lvert G_{ij}\cap T_k\rvert$ is computed with Shapely polygon clipping. An STRtree spatial
index limits the computation to candidate triangles.

### 13.2 Occupancy fraction

$$
\boxed{\phi_{ij}=\frac{\sum_k\lvert G_{ij}\cap T_k\rvert}{\lvert G_{ij}\rvert}}\in[0,1]
$$

$\phi_{ij}$ is the fraction of the cell covered by biofilm. Empty cells have $\phi=0$.

### 13.3 Area-weighted field averaging

For any element-averaged scalar $q_k$:

$$
\boxed{\overline q_{ij}=\frac{\sum_k\lvert G_{ij}\cap T_k\rvert\,q_k}{\sum_k\lvert G_{ij}\cap T_k\rvert}}
\qquad(\text{NaN where }\phi_{ij}=0).
$$

This is the average of $q$ over the *biofilm part* of the cell, not over the whole cell.

### 13.4 Conservation properties

Because the intersections partition each triangle exactly, provided the grid covers the whole mesh:

$$
\sum_{ij}\phi_{ij}\lvert G_{ij}\rvert=\sum_kA_k,\qquad
\sum_{ij}\phi_{ij}\lvert G_{ij}\rvert\,\overline q_{ij}=\sum_kA_k\,q_k=\int_\Omega q\,dA .
$$

Total biofilm area and area integrals of any mapped scalar are conserved to round-off. Point values, peaks and
gradients are **not** conserved, since averaging smooths them. Conserving area-integrated stress is not the
same as conserving force.

---

## 14. Summary of the computational methodology

1. **Geometry:** a seeded stochastic thickness $h(s)$ with fixed mean and RMS roughness, defined on the material
   arc length of an idealized rounded V-shaped bellows convolution. A hypothetical thinned
   ("regulated") twin is formed by uniform scaling, $h_r=0.5\,h_u$. This is a geometric scenario, not a model of
   gene regulation.
2. **Kinematics:** the bellows opens from $60^\circ$ to $60.5^\circ$ without stretching. Fold arc length is kept
   via $R_1=R_0(\pi-\theta_0)/(\pi-\theta_1)$, and material points are tracked by $s$, giving the base
   displacement $\mathbf u_b(s)$.
3. **Mesh:** a normal-offset curvilinear layer $\mathbf r(s,\eta)=\mathbf r_b+\eta h\mathbf n$, discretized into
   14 970 counter-clockwise P1 triangles on a uniform $s$-grid and checked for validity.
4. **Mechanics:** quasi-static, small-strain, homogeneous isotropic linear elasticity in plane strain, with
   illustrative $E=64$ Pa and $\nu=0.3$. The base is perfectly bonded and moves with the metal; all other
   surfaces are traction-free. Solved with SfePy in SI units.
5. **Indicators:** element-averaged von Mises stress, principal tensile stress (including $\sigma_{zz}$),
   principal tensile strain, and interfacial normal and shear traction. These are loading measures, not
   failure predictions.
6. **Experiments:** 20 paired morphologies, plus one-at-a-time sensitivity to $E$, $\nu$, $\Delta\theta$ and $R_0$.
7. **Integration:** a conservative, area-weighted, one-way projection of element fields onto a Cartesian grid.

---

## References

Literature (bibliographic details checked against Crossref):

- **[L1]** Picioreanu C, van Loosdrecht MCM, Heijnen JJ. Two-dimensional model of biofilm detachment caused by internal stress from liquid flow. *Biotechnology and Bioengineering* 72(2):205–218 (2001). https://doi.org/10.1002/1097-0290(20000120)72:2<205::AID-BIT9>3.0.CO;2-L
- **[L2]** Picioreanu C, Blauert F, Horn H, Wagner M. Determination of mechanical properties of biofilms by modelling the deformation measured using optical coherence tomography. *Water Research* 145:588–598 (2018). https://doi.org/10.1016/j.watres.2018.08.070
- **[L3]** Hohne DN, Younger JG, Solomon MJ. Flexible microfluidic device for mechanical property characterization of soft viscoelastic solids such as bacterial biofilms. *Langmuir* 25(13):7743–7751 (2009). https://doi.org/10.1021/la803413x
- **[L4]** Böl M, Ehret AE, Bolea Albero A, Hellriegel J, Krull R. Recent advances in mechanical characterisation of biofilm and their significance for material modelling. *Critical Reviews in Biotechnology* 33(2):145–171 (2013). https://doi.org/10.3109/07388551.2012.679250
- **[L9]** Stanley NR, Lazazzera BA. Defining the genetic differences between wild and domestic strains of *Bacillus subtilis* that affect poly-γ-DL-glutamic acid production and biofilm formation. *Molecular Microbiology* 57(4):1143–1158 (2005). https://doi.org/10.1111/j.1365-2958.2005.04746.x
- **[L10]** Yu Y, Yan F, Chen Y, Jin C, Guo J, Chai Y. Poly-γ-glutamic acids contribute to biofilm formation and plant root colonization in selected environmental isolates of *Bacillus subtilis*. *Frontiers in Microbiology* 7:1811 (2016). https://doi.org/10.3389/fmicb.2016.01811

Mission context (local copy in `../docs/lit/mech/`; bibliographic details **not independently verified**):

- **[R1]** Gilbert et al. (2020), ICES-2020-566, as cited in the team's modelling master plan
  (`biofilm_dynamics/iGEM_Integrated_Modelling_Implementation_Master_Plan.md`, §1.1). Describes biomass
  accumulation in bellows folds and its release associated with bellows cycling.

The full reference list, including viscoelasticity sources [L5–L8], is in
[VALIDATION_AND_LIMITATIONS.md](VALIDATION_AND_LIMITATIONS.md#references).
