"""
Structural FEM solver for irregular biofilm layers.

Physics:
    - 2D plane-strain linear elasticity
    - Homogeneous, isotropic material
    - Small-strain approximation
    - Quasi-static equilibrium

Supported loading:
    1. Prescribed surface traction
    2. Uniform substrate stretching/compression
    3. Combined traction and substrate deformation

Geometry:
    Biofilm occupies:
        0 <= x <= L
        0 <= y <= h(x)

Units:
    Input geometry: micrometres
    FEM coordinates: metres
    Young's modulus: Pa
    Traction: Pa
    Displacement: metres
    Stress: Pa
"""

import numpy as np

from sfepy.discrete.fem import Mesh, FEDomain, Field

from sfepy.discrete import (
    FieldVariable,
    Material,
    Integral,
    Equation,
    Equations,
    Problem,
    Function,
)

from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.terms import Term
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton

from .layer_mesh import build_layer_mesh


UM = 1e-6


# ============================================================
# SUBSTRATE DISPLACEMENT
# ============================================================

class SubstrateDisplacement:
    """
    Prescribed horizontal substrate displacement.

    Supported modes:
        uniform:
            ux = strain * (x - origin)

        sinusoidal:
            ux = amplitude * sin(2*pi*(x - origin)/wavelength)

    All coordinates and displacements are in metres.
    """

    def __init__(
        self,
        mode="uniform",
        strain=0.0,
        origin_m=0.0,
        amplitude_m=0.0,
        wavelength_m=None,
    ):
        self.mode = mode
        self.strain = float(strain)
        self.origin_m = float(origin_m)
        self.amplitude_m = float(amplitude_m)
        self.wavelength_m = wavelength_m

    def __call__(self, ts, coors, bc=None, problem=None):

        x = coors[:, 0]

        if self.mode == "uniform":
            return self.strain * (x - self.origin_m)

        elif self.mode == "sinusoidal":
            return self.amplitude_m * np.sin(
                2 * np.pi * (x - self.origin_m)
                / self.wavelength_m
            )

        raise ValueError(f"Unknown substrate mode: {self.mode}")


# ============================================================
# STRUCTURAL FEM SOLVER
# ============================================================

def solve_biofilm(
    x,
    h,
    ny=60,
    E=64.0,
    nu=0.3,
    traction_x=0.01,
    traction_y=0.0,
    substrate_strain=None,
    substrate_origin=None,
    substrate_mode="uniform",
    substrate_amplitude=None,
    substrate_wavelength=None,
):
    """
    Solve 2D plane-strain linear elasticity for an irregular biofilm.

    Parameters
    ----------
    x, h : ndarray
        Biofilm geometry in micrometres.

    ny : int
        Number of vertical mesh intervals.

    E : float
        Young's modulus (Pa).

    nu : float
        Poisson's ratio.

    traction_x, traction_y : float
        Prescribed surface traction components (Pa).

    substrate_strain : float or None
        Prescribed uniform horizontal substrate strain.

        None:
            Fixed substrate.

        Positive:
            Substrate stretching.

        Negative:
            Substrate compression.

    substrate_origin : float or None
        Reference x-coordinate in micrometres where substrate
        horizontal displacement equals zero.

        None:
            Use the midpoint of the substrate.

    Returns
    -------
    dict
        Mesh, displacement, strain, stress, and verification metrics.
    """

    # --------------------------------------------------------
    # 1. Validate input parameters
    # --------------------------------------------------------

    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)

    if E <= 0:
        raise ValueError("Young's modulus must be positive.")

    if not (-1.0 < nu < 0.5):
        raise ValueError(
            "Poisson's ratio must satisfy -1 < nu < 0.5."
        )

    if substrate_strain is not None:
        if not np.isfinite(substrate_strain):
            raise ValueError("Substrate strain must be finite.")

    # --------------------------------------------------------
    # 2. Generate triangular mesh
    # --------------------------------------------------------

    points_um, triangles, boundaries = build_layer_mesh(
        x,
        h,
        ny=ny,
    )

    points_m = points_um * UM

    mesh = Mesh.from_data(
        "biofilm_mesh",
        points_m,
        None,
        [triangles.astype(np.int32)],
        [np.zeros(len(triangles), dtype=np.int32)],
        ["2_3"],
    )

    # --------------------------------------------------------
    # 3. Define FEM domain and regions
    # --------------------------------------------------------

    domain = FEDomain(
        "biofilm_domain",
        mesh,
    )

    omega = domain.create_region(
        "Omega",
        "all",
    )

    # Bottom boundary: y = 0
    base = domain.create_region(
        "Base",
        "vertices in (y < 1e-10)",
        "facet",
    )

    # Irregular upper boundary: y = h(x)
    def select_top(coors, domain=None):

        expected_height = np.interp(
            coors[:, 0],
            x * UM,
            h * UM,
        )

        return np.where(
            np.abs(
                coors[:, 1] - expected_height
            ) < 1e-10
        )[0]

    top = domain.create_region(
        "Top",
        "vertices by select_top",
        "facet",
        functions={
            "select_top": select_top,
        },
    )

    # --------------------------------------------------------
    # 4. Define displacement field
    # --------------------------------------------------------

    field = Field.from_args(
        "displacement",
        np.float64,
        "vector",
        omega,
        approx_order=1,
    )

    # Unknown displacement
    u = FieldVariable(
        "u",
        "unknown",
        field,
    )

    # Test function
    v = FieldVariable(
        "v",
        "test",
        field,
        primary_var_name="u",
    )

    # --------------------------------------------------------
    # 5. Define material properties
    # --------------------------------------------------------

    # Plane-strain elastic stiffness matrix
    D = stiffness_from_youngpoisson(
        2,
        E,
        nu,
        plane="strain",
    )

    material = Material(
        "biofilm",
        D=D,
    )

    # --------------------------------------------------------
    # 6. Define substrate boundary conditions
    # --------------------------------------------------------

    if substrate_strain is None and substrate_mode == "uniform":
    
        fix_base = EssentialBC(
            "fix_base",
            base,
            {"u.all": 0.0},
        )
    
    else:
    
        if substrate_mode == "uniform":
    
            if substrate_strain is None:
                raise ValueError(
                    "Uniform deformation requires substrate_strain."
                )
    
            if substrate_origin is None:
                x0 = 0.5 * (np.min(x) + np.max(x)) * UM
            else:
                x0 = float(substrate_origin) * UM
    
            substrate_callback = SubstrateDisplacement(
                mode="uniform",
                strain=substrate_strain,
                origin_m=x0,
            )
    
        elif substrate_mode == "sinusoidal":
    
            if substrate_amplitude is None:
                raise ValueError("Specify substrate_amplitude.")
    
            if substrate_wavelength is None:
                raise ValueError("Specify substrate_wavelength.")
    
            if substrate_wavelength <= 0:
                raise ValueError(
                    "Substrate wavelength must be positive."
                )
    
            if substrate_origin is None:
                x0 = np.min(x) * UM
            else:
                x0 = float(substrate_origin) * UM
    
            substrate_callback = SubstrateDisplacement(
                mode="sinusoidal",
                origin_m=x0,
                amplitude_m=substrate_amplitude * UM,
                wavelength_m=substrate_wavelength * UM,
            )
    
        else:
            raise ValueError(
                f"Unknown substrate mode: {substrate_mode}"
            )
    
        displacement_function = Function(
            "substrate_displacement",
            substrate_callback,
        )
    
        fix_base = EssentialBC(
            "deforming_base",
            base,
            {
                "u.0": displacement_function,
                "u.1": 0.0,
            },
        )

    # --------------------------------------------------------
    # 7. Define surface traction
    # --------------------------------------------------------

    traction = Material(
        "traction",
        val=np.array([
            [traction_x],
            [traction_y],
        ]),
    )

    # --------------------------------------------------------
    # 8. Construct equilibrium equations
    # --------------------------------------------------------

    integral = Integral(
        "i",
        order=2,
    )

    # Internal elastic forces
    term_elastic = Term.new(
        "dw_lin_elastic(biofilm.D, v, u)",
        integral,
        omega,
        biofilm=material,
        v=v,
        u=u,
    )

    # External surface traction
    term_traction = Term.new(
        "dw_surface_ltr(traction.val, v)",
        integral,
        top,
        traction=traction,
        v=v,
    )

    # Internal forces = external forces
    equation = Equation(
        "elastic_equilibrium",
        term_elastic - term_traction,
    )

    equations = Equations([
        equation,
    ])

    # --------------------------------------------------------
    # 9. Create and solve FEM problem
    # --------------------------------------------------------

    problem = Problem(
        "biofilm_elasticity",
        equations=equations,
    )

    # Apply substrate boundary condition
    problem.time_update(
        ebcs=Conditions([
            fix_base,
        ])
    )

    # Linear solver
    ls = ScipyDirect({})

    # Newton interface
    nls = Newton(
        {
            "i_max": 1,
            "eps_a": 1e-10,
        },
        lin_solver=ls,
    )

    problem.set_solver(nls)

    # Solve equilibrium equations
    state = problem.solve(
        save_results=False,
    )

    # --------------------------------------------------------
    # 10. Extract displacement
    # --------------------------------------------------------

    displacement = state["u"]().reshape(
        -1,
        2,
    )

    displacement_magnitude = np.linalg.norm(
        displacement,
        axis=1,
    )

    # --------------------------------------------------------
    # 11. Calculate strain
    # --------------------------------------------------------

    strain = problem.evaluate(
        "ev_cauchy_strain.i.Omega(u)",
        mode="el_avg",
        integrals={
            "i": integral,
        },
    )

    strain_elements = strain.reshape(
        -1,
        3,
    )

    exx = strain_elements[:, 0]
    eyy = strain_elements[:, 1]
    exy = strain_elements[:, 2]

    # Principal strains
    mean_strain = (
        exx + eyy
    ) / 2.0

    radius = np.sqrt(
        ((exx - eyy) / 2.0)**2
        + exy**2
    )

    eps1 = mean_strain + radius
    eps2 = mean_strain - radius

    max_principal_strain = np.max(
        np.maximum(
            np.abs(eps1),
            np.abs(eps2),
        )
    )

    # --------------------------------------------------------
    # 12. Calculate stress
    # --------------------------------------------------------

    stress = problem.evaluate(
        "ev_cauchy_stress.i.Omega(biofilm.D, u)",
        mode="el_avg",
        integrals={
            "i": integral,
        },
        biofilm=material,
    )

    stress_elements = stress.reshape(
        -1,
        3,
    )

    sxx = stress_elements[:, 0]
    syy = stress_elements[:, 1]
    txy = stress_elements[:, 2]

    # Out-of-plane stress under plane strain
    szz = nu * (
        sxx + syy
    )

    # Von Mises equivalent stress
    sigma_e = np.sqrt(
        0.5 * (
            (sxx - syy)**2
            + (syy - szz)**2
            + (szz - sxx)**2
        )
        + 3.0 * txy**2
    )

    # --------------------------------------------------------
    # 13. Calculate element areas
    # --------------------------------------------------------

    P = points_m[triangles]

    AB = P[:, 1] - P[:, 0]
    AC = P[:, 2] - P[:, 0]

    areas = 0.5 * np.abs(
        AB[:, 0] * AC[:, 1]
        - AB[:, 1] * AC[:, 0]
    )

    # --------------------------------------------------------
    # 14. Calculate triangle centroids
    # --------------------------------------------------------

    centroids_um = points_um[
        triangles
    ].mean(axis=1)

    xc = centroids_um[:, 0]
    yc = centroids_um[:, 1]

    # --------------------------------------------------------
    # 15. Calculate interior stress
    # --------------------------------------------------------

    # Exclude regions within 5 μm of bottom corners
    corner_radius = 5.0

    dist_left = np.hypot(
        xc - np.min(x),
        yc,
    )

    dist_right = np.hypot(
        xc - np.max(x),
        yc,
    )

    interior_mask = (
        (dist_left > corner_radius)
        & (dist_right > corner_radius)
    )

    global_peak = float(
        sigma_e.max()
    )

    interior_peak = float(
        sigma_e[interior_mask].max()
    )

    # --------------------------------------------------------
    # 16. Calculate applied surface force
    # --------------------------------------------------------

    top_points = points_m[
        boundaries["top"]
    ]

    edge_vectors = np.diff(
        top_points,
        axis=0,
    )

    surface_length = np.linalg.norm(
        edge_vectors,
        axis=1,
    ).sum()

    F_applied_x = (
        traction_x * surface_length
    )

    F_applied_y = (
        traction_y * surface_length
    )

    # --------------------------------------------------------
    # 17. Return simulation results
    # --------------------------------------------------------

    return {

        # Geometry
        "points_um": points_um,
        "points_m": points_m,
        "triangles": triangles,
        "boundaries": boundaries,

        # Displacement
        "displacement": displacement,
        "displacement_magnitude": displacement_magnitude,

        # Strain
        "strain": strain_elements,
        "max_principal_strain": float(
            max_principal_strain
        ),

        # Stress
        "stress": stress_elements,
        "sigma_e": sigma_e,
        "global_peak": global_peak,
        "interior_peak": interior_peak,

        # Geometry metrics
        "areas": areas,
        "centroids_um": centroids_um,
        "surface_length_um": float(
            surface_length / UM
        ),

        # Displacement metrics
        "max_displacement_um": float(
            displacement_magnitude.max() / UM
        ),

        # External force
        "F_applied_x": float(F_applied_x),
        "F_applied_y": float(F_applied_y),

        # Simulation parameters
        "E": float(E),
        "nu": float(nu),
        "traction_x": float(traction_x),
        "traction_y": float(traction_y),
        "substrate_strain": substrate_strain,
        "substrate_origin": substrate_origin,
    }