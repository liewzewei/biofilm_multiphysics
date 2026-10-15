"""
Structural FEM for biofilm attached to a curved bellows surface.

Physics:
    - 2D plane-strain linear elasticity
    - Quasi-static equilibrium
    - Prescribed substrate displacement
    - Traction-free exposed boundaries
    - Perfect biofilm-substrate adhesion

Units:
    Input coordinates: micrometres
    Input displacements: micrometres
    FEM coordinates: metres
    Young's modulus: Pa
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


UM = 1e-6


class PrescribedBellowsDisplacement:
    """
    Interpolate prescribed displacement along the curved substrate.

    The callback receives coordinates in metres and returns
    displacement values in metres.
    """

    def __init__(self, substrate_points, displacement):

        self.points = np.asarray(
            substrate_points,
            dtype=float,
        )

        self.displacement = np.asarray(
            displacement,
            dtype=float,
        )

        if self.points.shape != self.displacement.shape:
            raise ValueError(
                "Substrate points and displacements must have "
                "matching shapes."
            )

        if self.points.ndim != 2 or self.points.shape[1] != 2:
            raise ValueError(
                "Expected arrays of shape (N, 2)."
            )

    def _lookup(self, coors):
        """
        Match boundary coordinates to prescribed substrate nodes.
        """

        from scipy.spatial import cKDTree

        tree = cKDTree(self.points)

        distances, indices = tree.query(coors)

        if np.max(distances) > 1e-10:
            raise ValueError(
                "Boundary coordinates do not match substrate nodes."
            )

        return self.displacement[indices]

    def horizontal(self, ts, coors, bc=None, problem=None):
        return self._lookup(coors)[:, 0]

    def vertical(self, ts, coors, bc=None, problem=None):
        return self._lookup(coors)[:, 1]


def solve_curved_biofilm(
    points,
    triangles,
    boundaries,
    substrate_displacement,
    E=64.0,
    nu=0.3,
):
    """
    Solve plane-strain elasticity for a curved biofilm.

    Parameters
    ----------
    points : ndarray (N, 2)
        Mesh coordinates in micrometres.

    triangles : ndarray (M, 3)
        Triangular element connectivity.

    boundaries : dict
        Must contain "substrate" node indices.

    substrate_displacement : ndarray (N_base, 2)
        Prescribed substrate displacement in micrometres.

    E : float
        Young's modulus (Pa).

    nu : float
        Poisson's ratio.

    Returns
    -------
    dict
        Displacement, strain, stress and mesh information.
    """

    # --------------------------------------------------
    # 1. Validate inputs
    # --------------------------------------------------

    points = np.asarray(points, dtype=float)

    triangles = np.asarray(
        triangles,
        dtype=np.int32,
    )

    base_nodes = np.asarray(
        boundaries["substrate"],
        dtype=np.int32,
    )

    substrate_displacement = np.asarray(
        substrate_displacement,
        dtype=float,
    )

    if substrate_displacement.shape != (len(base_nodes), 2):
        raise ValueError(
            "Substrate displacement must have shape (N_base, 2)."
        )

    if E <= 0:
        raise ValueError("Young's modulus must be positive.")

    if not (-1.0 < nu < 0.5):
        raise ValueError("Invalid Poisson's ratio.")

    # --------------------------------------------------
    # 2. Convert coordinates to metres
    # --------------------------------------------------

    points_m = points * UM

    displacement_m = substrate_displacement * UM

    # --------------------------------------------------
    # 3. Construct SfePy mesh
    # --------------------------------------------------

    mesh = Mesh.from_data(
        "curved_biofilm",
        points_m,
        None,
        [triangles],
        [np.zeros(len(triangles), dtype=np.int32)],
        ["2_3"],
    )

    domain = FEDomain(
        "curved_domain",
        mesh,
    )

    omega = domain.create_region(
        "Omega",
        "all",
    )

    # --------------------------------------------------
    # 4. Identify curved substrate
    # --------------------------------------------------

    substrate_points_m = points_m[base_nodes]

    def select_substrate(coors, domain=None):

        from scipy.spatial import cKDTree

        tree = cKDTree(substrate_points_m)

        distances, _ = tree.query(coors)

        return np.where(
            distances < 1e-10
        )[0]

    substrate = domain.create_region(
        "Substrate",
        "vertices by select_substrate",
        "facet",
        functions={
            "select_substrate": select_substrate,
        },
    )

    # --------------------------------------------------
    # 5. Define displacement field
    # --------------------------------------------------

    field = Field.from_args(
        "displacement",
        np.float64,
        "vector",
        omega,
        approx_order=1,
    )

    u = FieldVariable(
        "u",
        "unknown",
        field,
    )

    v = FieldVariable(
        "v",
        "test",
        field,
        primary_var_name="u",
    )

    # --------------------------------------------------
    # 6. Define elastic material
    # --------------------------------------------------

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

    # --------------------------------------------------
    # 7. Define prescribed displacement
    # --------------------------------------------------

    callback = PrescribedBellowsDisplacement(
        substrate_points_m,
        displacement_m,
    )

    ux_function = Function(
        "bellows_ux",
        callback.horizontal,
    )

    uy_function = Function(
        "bellows_uy",
        callback.vertical,
    )

    deform_substrate = EssentialBC(
        "deform_substrate",
        substrate,
        {
            "u.0": ux_function,
            "u.1": uy_function,
        },
    )

    # --------------------------------------------------
    # 8. Construct elasticity equation
    # --------------------------------------------------

    integral = Integral(
        "i",
        order=2,
    )

    term_elastic = Term.new(
        "dw_lin_elastic(biofilm.D, v, u)",
        integral,
        omega,
        biofilm=material,
        v=v,
        u=u,
    )

    equation = Equation(
        "elastic_equilibrium",
        term_elastic,
    )

    equations = Equations([
        equation,
    ])

    # --------------------------------------------------
    # 9. Solve
    # --------------------------------------------------

    problem = Problem(
        "curved_biofilm_elasticity",
        equations=equations,
    )

    problem.time_update(
        ebcs=Conditions([
            deform_substrate,
        ])
    )

    ls = ScipyDirect({})

    nls = Newton(
        {
            "i_max": 1,
            "eps_a": 1e-10,
        },
        lin_solver=ls,
    )

    problem.set_solver(nls)

    state = problem.solve(
        save_results=False,
    )

    # --------------------------------------------------
    # 10. Extract displacement
    # --------------------------------------------------

    displacement = state["u"]().reshape(
        -1,
        2,
    )

    displacement_magnitude = np.linalg.norm(
        displacement,
        axis=1,
    )

    # --------------------------------------------------
    # 11. Calculate strain
    # --------------------------------------------------

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

    # SfePy stores engineering shear strain gamma_xy = 2 * eps_xy;
    # the principal-strain formula needs the tensor component eps_xy
    exy = 0.5 * strain_elements[:, 2]

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
    # --------------------------------------------------
    # Maximum principal tensile strain
    # --------------------------------------------------
    
    principal_tensile_strain = np.maximum(
        eps1,
        0.0,
    )
    
    max_tensile_strain = float(
        principal_tensile_strain.max()
    )

    # --------------------------------------------------
    # 12. Calculate stress
    # --------------------------------------------------

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

    # Plane-strain out-of-plane stress
    szz = nu * (
        sxx + syy
    )

    # --------------------------------------------------
    # Maximum principal tensile stress
    # --------------------------------------------------
    
    mean_stress = (sxx + syy) / 2.0
    
    stress_radius = np.sqrt(
        ((sxx - syy) / 2.0)**2
        + txy**2
    )
    
    # In-plane principal stresses
    sigma1 = mean_stress + stress_radius
    sigma2 = mean_stress - stress_radius
    
    # Maximum positive principal stress, including
    # the out-of-plane stress under plane strain
    principal_tensile_stress = np.maximum.reduce([
        sigma1,
        sigma2,
        szz,
        np.zeros_like(sigma1),
    ])
    
    max_tensile_stress = float(
        principal_tensile_stress.max()
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

    # --------------------------------------------------
    # 13. Calculate element areas
    # --------------------------------------------------

    P = points_m[triangles]

    AB = P[:, 1] - P[:, 0]
    AC = P[:, 2] - P[:, 0]

    areas = 0.5 * np.abs(
        AB[:, 0] * AC[:, 1]
        - AB[:, 1] * AC[:, 0]
    )

    # --------------------------------------------------
    # 14. Return results
    # --------------------------------------------------

    return {
        "points_um": points,
        "points_m": points_m,
        "triangles": triangles,
        "boundaries": boundaries,

        "displacement": displacement,
        "displacement_magnitude": displacement_magnitude,

        "strain": strain_elements,
        "stress": stress_elements,
        "sigma_e": sigma_e,

        "areas": areas,

        "max_displacement_um": float(
            displacement_magnitude.max() / UM
        ),

        "max_principal_strain": float(
            max_principal_strain
        ),

        "global_peak": float(
            sigma_e.max()
        ),

        "E": float(E),
        "nu": float(nu),

        # Principal tensile indicators
        "principal_tensile_stress": principal_tensile_stress,
        "principal_tensile_strain": principal_tensile_strain,
        
        "max_tensile_stress": max_tensile_stress,
        "max_tensile_strain": max_tensile_strain,
    }