import numpy as np

from src.layer_mesh import build_layer_mesh
from src.mesh_quality import evaluate_mesh_quality


def test_flat_mesh_quality():

    x = np.linspace(0, 100, 51)
    h = np.full_like(x, 40.0)

    points, triangles, _ = build_layer_mesh(
        x, h, ny=20
    )

    q = evaluate_mesh_quality(points, triangles)

    assert np.all(q["signed_area"] > 0)

    assert np.all(np.isfinite(q["min_angle"]))

    assert np.all(np.isfinite(q["aspect_ratio"]))

    assert q["min_angle"].min() > 0


def test_irregular_mesh_quality():

    x = np.linspace(0, 100, 51)

    h = 40 + 10 * np.sin(
        2 * np.pi * x / 100
    )

    points, triangles, _ = build_layer_mesh(
        x, h, ny=20
    )

    q = evaluate_mesh_quality(points, triangles)

    assert np.all(q["signed_area"] > 0)

    assert q["min_angle"].min() > 0