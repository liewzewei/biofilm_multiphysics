import numpy as np


def build_curved_layer_mesh(
    substrate_points,
    thickness,
    ny=40,
):
    """
    Construct a triangular biofilm mesh attached to a curved substrate.

    Parameters
    ----------
    substrate_points : ndarray (N, 2)
        Ordered bellows surface coordinates in micrometres.

    thickness : ndarray (N,)
        Biofilm thickness at each substrate point in micrometres.

    ny : int
        Number of intervals through biofilm thickness.

    Returns
    -------
    points : ndarray (N_nodes, 2)
        Mesh coordinates in micrometres.

    triangles : ndarray (N_triangles, 3)
        Counterclockwise triangle connectivity.

    boundaries : dict
        Node indices for substrate, exposed surface, and lateral ends.
    """

    substrate_points = np.asarray(
        substrate_points, dtype=float
    )

    thickness = np.asarray(thickness, dtype=float)

    n = len(substrate_points)

    if substrate_points.shape != (n, 2):
        raise ValueError("Expected substrate coordinates of shape (N, 2).")

    if thickness.shape != (n,):
        raise ValueError("Thickness must have one value per substrate point.")

    if np.any(thickness <= 0):
        raise ValueError("Biofilm thickness must be positive.")

    # --------------------------------------------------
    # 1. Calculate surface tangents
    # --------------------------------------------------

    dx = np.gradient(substrate_points[:, 0])
    dy = np.gradient(substrate_points[:, 1])

    tangent = np.column_stack([dx, dy])

    tangent /= np.linalg.norm(
        tangent, axis=1
    )[:, None]

    # --------------------------------------------------
    # 2. Calculate surface normals
    # --------------------------------------------------

    normals = np.column_stack([
        -tangent[:, 1],
        tangent[:, 0],
    ])

    # --------------------------------------------------
    # 3. Generate nodes through thickness
    # --------------------------------------------------

    eta = np.linspace(0, 1, ny + 1)

    grid = (
        substrate_points[None, :, :]
        + eta[:, None, None]
        * thickness[None, :, None]
        * normals[None, :, :]
    )

    points = grid.reshape(-1, 2)

    # --------------------------------------------------
    # 4. Generate triangular connectivity
    # --------------------------------------------------

    triangles = []

    for j in range(ny):
        for i in range(n - 1):

            n00 = j * n + i
            n10 = j * n + i + 1
            n01 = (j + 1) * n + i
            n11 = (j + 1) * n + i + 1

            triangles.append([
                n00, n10, n11
            ])

            triangles.append([
                n00, n11, n01
            ])

    triangles = np.asarray(
        triangles, dtype=np.int32
    )

    # Ensure counterclockwise orientation
    P = points[triangles]

    signed_area = 0.5 * (
        (P[:, 1, 0] - P[:, 0, 0])
        * (P[:, 2, 1] - P[:, 0, 1])
        -
        (P[:, 1, 1] - P[:, 0, 1])
        * (P[:, 2, 0] - P[:, 0, 0])
    )

    inverted = signed_area < 0

    triangles[inverted] = triangles[inverted][:, [0, 2, 1]]

    # --------------------------------------------------
    # 5. Identify boundaries
    # --------------------------------------------------

    boundaries = {
        "substrate": np.arange(n),
        "top": ny * n + np.arange(n),
        "left": np.arange(ny + 1) * n,
        "right": np.arange(ny + 1) * n + (n - 1),
    }

    return points, triangles, boundaries