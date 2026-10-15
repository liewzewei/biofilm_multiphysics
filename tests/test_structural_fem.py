import numpy as np

from src.stochastic_geometry import random_height_profile
from src.structural_fem import solve_biofilm


def test_structural_fem_baseline():

    x = np.linspace(0, 500, 501)

    h = random_height_profile(
        x,
        mean_height=80.0,
        roughness=15.0,
        n_features=10,
        min_height=5.0,
        rng=np.random.default_rng(0),
    )

    result = solve_biofilm(
        x=x,
        h=h,
        ny=60,
        E=64.0,
        nu=0.3,
        traction_x=0.01,
    )

    assert np.isclose(
        result["global_peak"],
        0.068521,
        rtol=0.01,
    )

    assert np.isclose(
        result["interior_peak"],
        0.039652,
        rtol=0.01,
    )

    assert result["max_principal_strain"] < 0.01

    assert np.all(np.isfinite(result["sigma_e"]))