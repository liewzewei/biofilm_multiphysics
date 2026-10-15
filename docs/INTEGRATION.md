# Integrating the Structural FEM with the Biofilm Dispersion Model

**Audience:** teammates who maintain the 2D Cartesian-grid dispersion model.
**Entry point:** `fem_to_grid()` in [`src/grid_mapping.py`](../src/grid_mapping.py).

> **The current integration is one-way: FEM → Cartesian grid.**
> Mechanical fields computed on our triangular mesh can be projected onto any Cartesian grid. Nothing flows
> back: dispersion-model fields do not yet affect the FEM geometry, stiffness or loading. No bidirectional
> coupling between the biological PDEs and structural mechanics has been implemented.

---

## 1. Purpose of integration

The structural model ([METHODOLOGY.md](METHODOLOGY.md)) computes where a biofilm attached to a flexing bellows
fold is mechanically loaded: von Mises stress, principal tensile stress and strain. The dispersion model
computes how cells, γ-PGA, nutrients and pressure evolve on a Cartesian grid. The interface puts the mechanical
fields onto the dispersion grid, so that the two can be plotted, compared, or used in future coupling
experiments (§17) without changing either model's equations.

## 2. Overview of the two computational models

| | Structural FEM (this repository) | Dispersion model (teammates) |
|---|---|---|
| Physics | Quasi-static plane-strain linear elasticity | Transport/reaction PDEs |
| Fields | Displacement (nodal); strain and stress (per element) | Motile cell density, non-motile cell density, γ-PGA concentration, nutrient concentration, effective pressure, Darcy velocity |
| Discretization | Unstructured P1 triangles on a curved strip | Structured Cartesian grid |
| Time | Static snapshot (one bellows movement) | Time-dependent |
| Geometry | Biofilm in one bellows crevice (about 235 × 187 μm) | Set by the dispersion model |

**Do not equate the two "pressures".** The dispersion model's *effective pressure* is a transport variable that
drives Darcy flow in the biofilm–fluid mixture. FEM stresses are Cauchy stresses in an elastic solid.
Interfacial *tractions* are forces per area on a boundary. These are different physical quantities: neither
should be fed into the other as if they were equivalent.

## 3. Coordinate systems and physical units

| Quantity | FEM output | Unit |
|---|---|---|
| Node coordinates `result["points_um"]` | Crevice frame (below) | **μm** |
| Displacement `result["displacement"]` | per node | **m** |
| Stress `result["stress"]`, `sigma_e`, `principal_tensile_stress` | per element | **Pa** |
| Strain `result["strain"]`, `principal_tensile_strain` | per element | dimensionless |
| Element areas `result["areas"]` | per element | **m²** |

**FEM crevice frame.** $x$ points right and $y$ up. The bellows is symmetric about $x=0$, the fold apex (metal
surface) is at $y=0$, and the crevice opens upward. For the baseline geometry, the biofilm occupies roughly
$x\in[-117,117]$ μm and $y\in[0,187]$ μm. This frame is a modelling convenience and has no fixed relation to
the dispersion domain (see §16).

**Units for `fem_to_grid`.** The mapping is pure geometry. `points`, `x_edges` and `y_edges` must use the
**same length unit**; the documented convention is μm. Field values pass through in whatever unit they carry
(Pa for stress). If the dispersion grid is stored in metres, divide its edges by `UM = 1e-6`, or equivalently
multiply by $10^6$.

## 4. FEM mesh representation

- `points`: `(N_nodes, 2)` float array of node coordinates.
- `triangles`: `(M, 3)` integer array of node indices for each counter-clockwise triangle.
- Element fields: `(M,)` arrays with **one constant value per triangle**. This is exact for P1 strain and
  stress, which are piecewise constant.

Production mesh: 8 000 nodes and 14 970 triangles per biofilm.

## 5. Cartesian grid representation

The grid is described only by its **cell edges**:

- `x_edges`: `(nx + 1,)`, strictly increasing.
- `y_edges`: `(ny + 1,)`, strictly increasing.

Cell $(j,i)$ is the rectangle $[x_i,x_{i+1}]\times[y_j,y_{j+1}]$. Spacing does **not** have to be uniform. This
matches the team grid contract (`dx_m`, `dy_m`, `origin_m`, `axis_order="yx"`), with
`x_edges = origin_m[0] + dx_m * arange(nx + 1)` and similarly for $y$.

## 6. Function API and parameters

```python
from src.grid_mapping import fem_to_grid

grid = fem_to_grid(points, triangles, fields, x_edges, y_edges)
```

| Parameter | Type | Description |
|---|---|---|
| `points` | `ndarray (N, 2)` | FEM node coordinates, same length unit as the edges (μm by convention) |
| `triangles` | `ndarray (M, 3)` int | Triangle connectivity |
| `fields` | `dict[str, ndarray (M,)]` | Any number of **element-averaged scalar** fields, keyed by name |
| `x_edges`, `y_edges` | `ndarray` | Cell edges, strictly increasing, at least 2 entries each |

The function raises `ValueError` if the grid has no cells, if the edges are not strictly increasing, or if a
field does not have exactly one value per triangle.

**Fields suitable for mapping** (keys of the dict returned by `solve_curved_biofilm`):

| Result key | Meaning |
|---|---|
| `sigma_e` | Von Mises equivalent stress (Pa) |
| `principal_tensile_stress` | Maximum principal tensile stress (Pa) |
| `principal_tensile_strain` | Maximum principal tensile strain (dimensionless) |

Individual tensor components, such as `result["stress"][:, k]`, can also be mapped, but they are tied to the FEM
axes (§16).

## 7. Required input array shapes

| Array | Shape | dtype |
|---|---|---|
| `points` | `(N, 2)` | float |
| `triangles` | `(M, 3)` | int (cast to `int32` internally) |
| each `fields[name]` | `(M,)` | float |
| `x_edges` | `(nx + 1,)` | float |
| `y_edges` | `(ny + 1,)` | float |

Nodal fields, such as displacement, must first be converted to element values (for example by averaging the
three nodes of each triangle) before they can be mapped. This conversion is not provided.

## 8. Output dictionary structure

```python
{
    "occupancy_fraction": ndarray (ny, nx),          # φ in [0, 1]
    "fields": {name: ndarray (ny, nx), ...},         # area-weighted averages, NaN where φ = 0
    "x_edges": ndarray (nx + 1,),                    # copies of the inputs
    "y_edges": ndarray (ny + 1,),
}
```

## 9. Occupancy fraction

$$
\phi_{ij}=\frac{\sum_k\lvert G_{ij}\cap T_k\rvert}{\lvert G_{ij}\rvert}
$$

$\phi_{ij}$ is the fraction of cell $G_{ij}$ covered by biofilm triangles $T_k$, using exact Shapely polygon
intersections. Cells fully inside the biofilm have $\phi=1$, boundary cells $0<\phi<1$, and empty cells
$\phi=0$.

$\phi$ is a **geometric** occupancy of the FEM biofilm. It is not a cell density and is not derived from any
dispersion-model variable.

## 10. Area-weighted stress and strain mapping

$$
\overline q_{ij}=\frac{\sum_k\lvert G_{ij}\cap T_k\rvert\,q_k}{\sum_k\lvert G_{ij}\cap T_k\rvert}
$$

$\overline q_{ij}$ is the average of $q$ over the **biofilm-covered part** of the cell. To obtain an amount per
whole cell, for example the area integral in the cell, use $\phi_{ij}\lvert G_{ij}\rvert\,\overline q_{ij}$.
Averaging smooths peaks, so $\max_{ij}\overline q_{ij}\le\max_k q_k$.

## 11. Treatment of empty grid cells (NaN)

Where $\phi_{ij}=0$, every mapped field is `NaN` and occupancy is `0.0`. `NaN` means "no biofilm here", not
"zero stress". Recommended handling:

- Integrals and sums: `np.nansum(field * phi * cell_area)`.
- Plots: leave the NaN values as they are; matplotlib renders them transparent.
- A consumer that cannot accept NaN: `np.where(phi > 0, field, fill)`, choosing `fill` explicitly, and record
  that choice.

## 12. Grid orientation and indexing conventions

- Arrays are indexed `[j, i]`, i.e. `[row, column]` = `[y, x]`. This is `axis_order = "yx"`.
- `j = 0` is the row with the **smallest** $y$, the bottom of the grid. `i = 0` is the column with the smallest
  $x$.
- With `imshow`, use `origin="lower"` and `extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]` to plot in
  physical orientation.

## 13. Example: mapping onto a dispersion-model grid

This script was run end to end with the repository's `src/` modules (conda env `biofilm-fem`). It reported a
relative area error of `0.00e+00` and an integral error of `2.19e-16`. Run it from the repository root.

```python
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(".").resolve()          # biofilm_multiphysics/
sys.path.insert(0, str(PROJECT_ROOT))

from sfepy.base.base import output
output.set_output(quiet=True)

from src.stochastic_geometry import random_height_profile
from src.rounded_bellows import deform_rounded_bellows, resample_bellows
from src.curved_layer_mesh import build_curved_layer_mesh
from src.curved_structural_fem import solve_curved_biofilm, UM
from src.grid_mapping import fem_to_grid

# ---------------------------------------------------------------
# 1. Run one structural FEM case (baseline parameters, seed 42)
# ---------------------------------------------------------------
deformation = deform_rounded_bellows(
    segment_length=200.0, initial_angle_deg=60.0,
    final_angle_deg=60.5, initial_radius=20.0,
)
s_um = np.linspace(0.0, deformation["s_initial"][-1], 500)
ref_um = resample_bellows(deformation["initial_points"], deformation["s_initial"], s_um)
cur_um = resample_bellows(deformation["final_points"], deformation["s_final"], s_um)
level = ref_um[:, 1].max()                  # reflect into upward-opening crevice
ref_um[:, 1] = level - ref_um[:, 1]
cur_um[:, 1] = level - cur_um[:, 1]
shift = min(ref_um[:, 1].min(), cur_um[:, 1].min())
ref_um[:, 1] -= shift
cur_um[:, 1] -= shift

h_um = random_height_profile(s_um, mean_height=10.0, roughness=3.0,
                             n_features=10, min_height=2.0,
                             rng=np.random.default_rng(42))
points_um, triangles, boundaries = build_curved_layer_mesh(ref_um, h_um, ny=15, normal_side="left")
result = solve_curved_biofilm(points_um, triangles, boundaries,
                              substrate_displacement=cur_um - ref_um, E=64.0, nu=0.3)

# ---------------------------------------------------------------
# 2. Describe the dispersion-model grid (metres, axis order "yx")
#    Replace these four values with the dispersion model's own.
# ---------------------------------------------------------------
origin_m = np.array([-150e-6, -10e-6])      # lower-left corner (x0, y0)
dx_m, dy_m = 2e-6, 2e-6                     # cell size
nx, ny = 150, 110                           # cell counts

x_edges_m = origin_m[0] + dx_m * np.arange(nx + 1)
y_edges_m = origin_m[1] + dy_m * np.arange(ny + 1)

# ---------------------------------------------------------------
# 3. Place the FEM geometry in the grid's frame (rigid transform)
#    Here: identity rotation, zero offset. Edit to align frames.
# ---------------------------------------------------------------
theta = 0.0                                 # rotation (rad), counter-clockwise
offset_um = np.array([0.0, 0.0])            # translation (μm)
rot = np.array([[np.cos(theta), -np.sin(theta)],
                [np.sin(theta),  np.cos(theta)]])
points_grid_um = points_um @ rot.T + offset_um

# ---------------------------------------------------------------
# 4. Map element-averaged scalar fields (all inputs in μm)
# ---------------------------------------------------------------
grid = fem_to_grid(
    points=points_grid_um,
    triangles=result["triangles"],
    fields={
        "von_mises_Pa": result["sigma_e"],
        "tensile_stress_Pa": result["principal_tensile_stress"],
    },
    x_edges=x_edges_m / UM,                 # metres -> μm
    y_edges=y_edges_m / UM,
)

phi = grid["occupancy_fraction"]            # shape (ny, nx), j=0 is the bottom row
sigma_t = grid["fields"]["tensile_stress_Pa"]   # NaN where phi == 0

# ---------------------------------------------------------------
# 5. Conservation check (only valid if the grid covers the mesh)
# ---------------------------------------------------------------
cell_area_m2 = dx_m * dy_m
fem_area_m2 = result["areas"].sum()         # solver areas are in m²
grid_area_m2 = np.sum(phi) * cell_area_m2
I_fem = np.sum(result["principal_tensile_stress"] * result["areas"])
I_grid = np.nansum(sigma_t * phi) * cell_area_m2
print(f"area error:     {abs(grid_area_m2 - fem_area_m2) / fem_area_m2:.2e}")
print(f"integral error: {abs(I_grid - I_fem) / abs(I_fem):.2e}")

# ---------------------------------------------------------------
# 6. Hand-off arrays (NaN -> 0 only where a consumer requires it)
# ---------------------------------------------------------------
np.savez_compressed(
    "fem_on_dispersion_grid.npz",
    occupancy_fraction=phi,
    tensile_stress_Pa=sigma_t,
    von_mises_Pa=grid["fields"]["von_mises_Pa"],
    dx_m=dx_m, dy_m=dy_m, origin_m=origin_m, axis_order="yx",
)
```

The FEM part takes under a second, and mapping a 150 × 110 grid takes well under a second. The final notebook
(Section 16) shows the same mapping on a 100 × 100 grid with plots.

## 14. Conservation verification

When the grid covers the entire mesh, the mapping exactly preserves:

$$
\sum_{ij}\phi_{ij}\lvert G_{ij}\rvert=\sum_kA_k,
\qquad
\sum_{ij}\phi_{ij}\lvert G_{ij}\rvert\,\overline q_{ij}=\sum_kA_k\,q_k .
$$

Verified results: area error $0$, $\int\sigma_t\,dA$ error $2.4\times10^{-16}$, $\int\sigma_e\,dA$ error $0$
(final notebook, 100 × 100 grid), and $0$ / $2.2\times10^{-16}$ for the example above.

**Run step 5 of the example for every new grid.** An area error well above $10^{-12}$ means part of the mesh lies
outside the grid, or the units of the points and the edges do not match.

Remember: `result["areas"]` is in **m²**. Convert consistently, either both sides in m² as in the example, or
both in μm² by dividing the FEM areas by `UM**2`.

## 15. How to change grid resolution

Only `x_edges` and `y_edges` change. The FEM mesh and the mapping code are **independent of the dispersion grid
resolution**: no remeshing or re-solving is needed.

```python
dx_m = dy_m = 1e-6                        # finer grid
nx = int(np.ceil((x_max_m - origin_m[0]) / dx_m))
ny = int(np.ceil((y_max_m - origin_m[1]) / dy_m))
x_edges_m = origin_m[0] + dx_m * np.arange(nx + 1)
y_edges_m = origin_m[1] + dy_m * np.arange(ny + 1)
```

Guidance:

- **Cell size relative to the FEM resolution.** Cells much smaller than the FEM elements (about 0.9 μm along the
  surface, 0.2–1 μm through the thickness) add no information, because each element is constant. Cells much
  larger than the biofilm thickness (5–10 μm) smear the film into partially occupied cells and flatten stress
  peaks.
- **Cost.** The mapping loops over all cells in Python. Cost grows with `nx * ny`, plus the number of triangles
  that overlap each cell. If the dispersion domain is much larger than the crevice, map onto the sub-grid that
  covers the biofilm's bounding box, then insert the result into the full arrays (occupancy 0, fields NaN
  elsewhere).

## 16. How to align our coordinate system with the dispersion model

The FEM frame (§3) is arbitrary up to a rigid motion. To align it:

1. **Choose a correspondence.** Decide where the bellows surface lies in the dispersion domain. For example,
   place the fold apex at a specific dispersion-grid point and align the crevice axis with the grid's $y$ axis.
2. **Apply a rigid transform to the FEM node coordinates only:**
   $$
   \mathbf x_{\text{grid}}=\mathbf R(\vartheta)\,\mathbf x_{\text{FEM}}+\mathbf b,
   \qquad
   \mathbf R(\vartheta)=\begin{bmatrix}\cos\vartheta&-\sin\vartheta\\ \sin\vartheta&\cos\vartheta\end{bmatrix}.
   $$
   In the example these are `theta` and `offset_um`. Rigid transforms preserve areas, so conservation still
   holds.
3. **Scalar invariants are unaffected.** $\sigma_e$, $\sigma_t$, $\varepsilon_t$ and $\phi$ do not depend on the
   choice of axes, so they can be mapped after any rotation.
4. **Tensor components are affected.** $\sigma_{xx}$, $\sigma_{yy}$ and $\tau_{xy}$ are defined in FEM axes. After
   a rotation by $\vartheta$, transform them with $\boldsymbol\sigma'=\mathbf R\boldsymbol\sigma\mathbf R^{T}$
   before mapping, or map only invariants.
5. **Do not rescale the geometry** to fit the grid. Changing lengths would change the physics; refine the grid
   instead.
6. **Mirror images.** If the dispersion domain's "substrate side" is at the top, use a reflection
   ($y\rightarrow-y$) instead of a rotation. The triangle orientation flips under reflection, but `fem_to_grid` uses
   unsigned intersection areas, so the mapping is unaffected.

Record `theta`, `offset_um` and any reflection with the exported arrays.

## 17. Potential future coupling (not implemented)

The coupling paths below are **proposals**. Each requires constitutive information we do not yet have.

| Direction | Idea | What would be needed |
|---|---|---|
| Biomass density → geometry | Extract the biofilm outline from a thresholded biomass field and build the FEM mesh from it, replacing the synthetic $h(s)$ | A justified occupancy threshold; an inverse mapping from grid to boundary; handling of holes, islands and overhangs, which the current normal-offset mesher cannot represent |
| γ-PGA concentration → stiffness | Spatially varying $E(\mathbf x)=f(C_{\gamma\text{PGA}})$ | Measured modulus–concentration data. The FEM currently supports only uniform $E$, so `solve_curved_biofilm` would need per-element material values |
| γ-PGA concentration → strength | Local strength $\sigma_c(\mathbf x)=g(C_{\gamma\text{PGA}})$ used to normalize stress ($\sigma_t/\sigma_c$) | Measured cohesive or tensile strength; a validated failure criterion |
| Mechanical loading → fragmentation | Use mapped $\sigma_t$ or interface $t_n^+$ as an input to a detachment rate in the dispersion model | A calibrated detachment law. **Stress alone does not predict fragmentation, fragment size or detached biomass** |

Even with these links, the dispersion model's effective pressure should not be applied as a mechanical
traction without a derivation that connects the two (§2).

## 18. Limitations of the current integration

1. **One-way, static.** The mapping transfers a single FEM snapshot (one bellows movement). It has no time
   dependence and does not feed back.
2. **Geometry is not shared.** The FEM biofilm comes from the synthetic generator, not from the dispersion
   model's biomass field. Overlaying them on one grid does not mean they describe the same biofilm.
3. **Scalar element fields only.** No nodal fields, no interface tractions (boundary quantities), and no
   automatic tensor rotation.
4. **Smoothing.** Cell averages under-represent peaks. Use FEM values directly for peak-stress statements.
5. **Coverage.** Conservation holds only if the grid covers the whole mesh. The function does not warn when
   part of the mesh falls outside the grid.
6. **Inherited model limitations.** Everything in
   [VALIDATION_AND_LIMITATIONS.md §10](VALIDATION_AND_LIMITATIONS.md#10-scientific-limitations) applies to the
   mapped fields, including illustrative $E$ and $\nu$, and the absence of viscoelasticity and fracture.
