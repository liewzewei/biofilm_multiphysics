import numpy as np


def build_layer_mesh(x, h, ny=20):
    """
    Generate a boundary-fitted triangular mesh for a biofilm layer.

    Parameters
    ----------
    x : 1D array
        Horizontal coordinates.
    h : 1D array
        Biofilm height at each x-coordinate.
    ny : int
        Number of vertical element intervals.

    Returns
    -------
    points : ndarray, shape (n_nodes, 2)
        Physical coordinates of mesh nodes.
    triangles : ndarray, shape (n_triangles, 3)
        Node indices defining each triangle.
    boundaries : dict
        Node indices along the four boundaries.
    """

    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)

    # Validate geometry
    if x.ndim != 1 or h.ndim != 1:
        raise ValueError("x and h must be 1D arrays")

    if len(x) != len(h):
        raise ValueError("x and h must have the same length")

    if len(x) < 2 or ny < 1:
        raise ValueError("Need at least 2 x-points and 1 vertical interval")

    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(h)):
        raise ValueError("Coordinates must be finite")

    if not np.all(np.diff(x) > 0):
        raise ValueError("x must be strictly increasing")

    if not np.all(h > 0):
        raise ValueError("Biofilm height must be positive")

    # Number of horizontal nodes
    nx = len(x)

    # Normalized vertical coordinates
    eta = np.linspace(0, 1, ny + 1)

    # Construct physical coordinate arrays
    X = np.broadcast_to(x[None, :], (ny + 1, nx))
    Y = eta[:, None] * h[None, :]

    # Flatten into (x, y) coordinate pairs
    points = np.column_stack([
        X.ravel(),
        Y.ravel()
    ])

    # Convert 2D indices to flattened node IDs
    def node_id(i, j):
        return j * nx + i

    # Generate triangles
    triangles = []

    for j in range(ny):
        for i in range(nx - 1):

            n00 = node_id(i, j)
            n10 = node_id(i + 1, j)
            n01 = node_id(i, j + 1)
            n11 = node_id(i + 1, j + 1)

            # Alternate diagonals to reduce directional bias
            if (i + j) % 2 == 0:
                triangles.append([n00, n10, n11])
                triangles.append([n00, n11, n01])
            else:
                triangles.append([n00, n10, n01])
                triangles.append([n10, n11, n01])

    triangles = np.asarray(triangles, dtype=int)

    # Identify boundary nodes
    boundaries = {
        "base": np.array(
            [node_id(i, 0) for i in range(nx)]
        ),
        "top": np.array(
            [node_id(i, ny) for i in range(nx)]
        ),
        "left": np.array(
            [node_id(0, j) for j in range(ny + 1)]
        ),
        "right": np.array(
            [node_id(nx - 1, j) for j in range(ny + 1)]
        )
    }

    return points, triangles, boundaries