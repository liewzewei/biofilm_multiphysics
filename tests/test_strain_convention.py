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

from src.layer_mesh import build_layer_mesh
from src.curved_layer_mesh import build_curved_layer_mesh
from src.curved_structural_fem import solve_curved_biofilm
from src.structural_fem import solve_biofilm


def test_sfepy_returns_engineering_shear_strain():
    # Homogeneous pure shear u = (a y, a x) on the whole boundary
    # gives tensor shear strain eps_xy = a everywhere
    a = 1e-3

    x = np.linspace(0, 10, 11)
    points, triangles, _ = build_layer_mesh(x, np.full_like(x, 10.0), ny=10)

    mesh = Mesh.from_data(
        "square",
        points,
        None,
        [triangles.astype(np.int32)],
        [np.zeros(len(triangles), dtype=np.int32)],
        ["2_3"],
    )
    domain = FEDomain("square", mesh)
    omega = domain.create_region("Omega", "all")
    boundary = domain.create_region("Boundary", "vertices of surface", "facet")

    field = Field.from_args("u", np.float64, "vector", omega, approx_order=1)
    u = FieldVariable("u", "unknown", field)
    v = FieldVariable("v", "test", field, primary_var_name="u")

    material = Material("m", D=stiffness_from_youngpoisson(2, 64.0, 0.3, plane="strain"))
    integral = Integral("i", order=2)

    ux = Function("ux", lambda ts, coors, **kwargs: a * coors[:, 1])
    uy = Function("uy", lambda ts, coors, **kwargs: a * coors[:, 0])

    problem = Problem(
        "pure_shear",
        equations=Equations([
            Equation(
                "balance",
                Term.new("dw_lin_elastic(m.D, v, u)", integral, omega, m=material, v=v, u=u),
            )
        ]),
    )
    problem.time_update(
        ebcs=Conditions([EssentialBC("shear", boundary, {"u.0": ux, "u.1": uy})])
    )
    problem.set_solver(Newton({"i_max": 1, "eps_a": 1e-10}, lin_solver=ScipyDirect({})))
    problem.solve(save_results=False)

    strain = problem.evaluate(
        "ev_cauchy_strain.i.Omega(u)", mode="el_avg", integrals={"i": integral}
    ).reshape(-1, 3)

    # Third component is gamma_xy = 2 eps_xy, which the solvers must halve
    assert np.allclose(strain[:, 2], 2 * a, rtol=1e-10)


def _solve_strip(rotation_deg):
    # Biofilm strip on a straight substrate, loaded by a nonuniform
    # substrate displacement; the whole problem is rotated rigidly
    s = np.linspace(0.0, 100.0, 41)
    substrate = np.column_stack([s, np.zeros_like(s)])
    displacement = np.column_stack([
        0.5 * np.sin(2 * np.pi * s / 100.0),
        0.3 * np.cos(2 * np.pi * s / 50.0),
    ])

    theta = np.deg2rad(rotation_deg)
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)],
    ])

    points, triangles, boundaries = build_curved_layer_mesh(
        substrate @ R.T,
        np.full(len(s), 5.0),
        ny=4,
    )

    return solve_curved_biofilm(
        points,
        triangles,
        boundaries,
        substrate_displacement=displacement @ R.T,
    )


def test_curved_principal_strain_is_rotation_invariant():
    # Principal strains are tensor invariants: rotating the geometry and the
    # prescribed displacement together must not change them. Using
    # gamma_xy in place of eps_xy breaks this invariance.
    reference = _solve_strip(0.0)
    rotated = _solve_strip(30.0)

    assert reference["max_principal_strain"] > 1e-4

    assert np.allclose(
        rotated["principal_tensile_strain"],
        reference["principal_tensile_strain"],
        rtol=1e-8,
        atol=1e-14,
    )
    assert np.isclose(
        rotated["max_principal_strain"],
        reference["max_principal_strain"],
        rtol=1e-8,
    )

    # Stress invariants were already correct and must remain so
    assert np.allclose(rotated["sigma_e"], reference["sigma_e"], rtol=1e-8, atol=1e-14)


def test_flat_principal_strain_matches_tensor_eigenvalues():
    x = np.linspace(0, 100, 51)
    h = 40 + 10 * np.sin(2 * np.pi * x / 100)

    result = solve_biofilm(x=x, h=h, ny=10, traction_x=0.01)

    strain = result["strain"]

    # Tensor built from SfePy's engineering shear strain
    eps = np.zeros((len(strain), 2, 2))
    eps[:, 0, 0] = strain[:, 0]
    eps[:, 1, 1] = strain[:, 1]
    eps[:, 0, 1] = eps[:, 1, 0] = 0.5 * strain[:, 2]

    eigenvalues = np.linalg.eigvalsh(eps)

    assert np.isclose(
        result["max_principal_strain"],
        np.abs(eigenvalues).max(),
        rtol=1e-10,
    )
