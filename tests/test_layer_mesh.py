import numpy as np

from src.layer_mesh import build_layer_mesh


def test_flat_layer():
    x = np.linspace(0, 100, 51)
    h = np.full_like(x, 40.0)

    points, triangles, boundaries = build_layer_mesh(
        x, h, ny=20
    )

    # Expected number of nodes and triangles
    assert len(points) == 51 * 21
    assert len(triangles) == 50 * 20 * 2

    # Base must lie at y = 0
    assert np.allclose(
        points[boundaries["base"], 1],
        0
    )

    # Top must lie at y = 40
    assert np.allclose(
        points[boundaries["top"], 1],
        40
    )


def test_mesh_area():
    x = np.linspace(0, 100, 51)
    h = 40 + 10 * np.sin(2 * np.pi * x / 100)

    points, triangles, _ = build_layer_mesh(
        x, h, ny=20
    )

    # Coordinates of every triangle's vertices
    P = points[triangles]

    AB = P[:, 1] - P[:, 0]
    AC = P[:, 2] - P[:, 0]

    signed_areas = 0.5 * (
        AB[:, 0] * AC[:, 1]
        - AB[:, 1] * AC[:, 0]
    )

    # All triangles must be counterclockwise
    assert np.all(signed_areas > 0)

    # Mesh area must match area under height profile
    area_mesh = signed_areas.sum()
    area_profile = np.trapezoid(h, x)

    assert np.isclose(area_mesh, area_profile)