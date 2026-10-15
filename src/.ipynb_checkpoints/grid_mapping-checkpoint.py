"""
Map element-averaged FEM fields onto a Cartesian grid.

The mapping uses exact polygon intersection areas.

Inputs:
    FEM coordinates: micrometres
    Grid coordinates: micrometres
    Fields: element-averaged values

Outputs:
    Cartesian arrays with shape (ny, nx)

Empty cells:
    Occupancy = 0
    Field values = NaN
"""

import numpy as np

from shapely.geometry import Polygon, box
from shapely.strtree import STRtree


def fem_to_grid(
    points,
    triangles,
    fields,
    x_edges,
    y_edges,
):
    """
    Map triangular FEM fields to a Cartesian grid.

    Parameters
    ----------
    points : ndarray (N, 2)
        FEM node coordinates in micrometres.

    triangles : ndarray (M, 3)
        Triangle connectivity.

    fields : dict
        Dictionary of element-averaged fields.

        Each array must have length M.

    x_edges, y_edges : ndarray
        Cartesian grid cell boundaries in micrometres.

    Returns
    -------
    dict
        Occupancy fraction, mapped fields,
        and grid coordinate arrays.
    """

    points = np.asarray(points, dtype=float)

    triangles = np.asarray(
        triangles,
        dtype=np.int32,
    )

    x_edges = np.asarray(x_edges, dtype=float)
    y_edges = np.asarray(y_edges, dtype=float)

    nx = len(x_edges) - 1
    ny = len(y_edges) - 1

    if nx <= 0 or ny <= 0:
        raise ValueError("Grid must contain at least one cell.")

    if np.any(np.diff(x_edges) <= 0):
        raise ValueError("x_edges must be strictly increasing.")

    if np.any(np.diff(y_edges) <= 0):
        raise ValueError("y_edges must be strictly increasing.")

    # --------------------------------------------------
    # 1. Construct FEM triangle polygons
    # --------------------------------------------------

    triangle_polygons = [
        Polygon(points[tri])
        for tri in triangles
    ]

    # Spatial index
    tree = STRtree(triangle_polygons)

    # --------------------------------------------------
    # 2. Validate fields
    # --------------------------------------------------

    fields = {
        name: np.asarray(values, dtype=float)
        for name, values in fields.items()
    }

    for name, values in fields.items():

        if values.shape != (len(triangles),):
            raise ValueError(
                f"Field '{name}' must have one value per triangle."
            )

    # --------------------------------------------------
    # 3. Initialize output arrays
    # --------------------------------------------------

    occupancy = np.zeros(
        (ny, nx),
        dtype=float,
    )

    mapped_fields = {
        name: np.full(
            (ny, nx),
            np.nan,
            dtype=float,
        )
        for name in fields
    }

    # --------------------------------------------------
    # 4. Calculate cell intersections
    # --------------------------------------------------

    for j in range(ny):

        for i in range(nx):

            cell = box(
                x_edges[i],
                y_edges[j],
                x_edges[i + 1],
                y_edges[j + 1],
            )

            cell_area = cell.area

            # Shapely 2.x returns integer indices
            candidate_indices = tree.query(cell)

            if len(candidate_indices) == 0:
                continue

            total_overlap = 0.0

            weighted_sums = {
                name: 0.0
                for name in fields
            }

            for k in candidate_indices:

                k = int(k)

                overlap = cell.intersection(
                    triangle_polygons[k]
                ).area

                if overlap <= 0:
                    continue

                total_overlap += overlap

                for name, values in fields.items():

                    weighted_sums[name] += (
                        overlap * values[k]
                    )

            if total_overlap <= 0:
                continue

            # Fraction of cell occupied by biofilm
            occupancy[j, i] = (
                total_overlap / cell_area
            )

            # Area-weighted field averages
            for name in fields:

                mapped_fields[name][j, i] = (
                    weighted_sums[name] / total_overlap
                )

    # --------------------------------------------------
    # 5. Return mapped results
    # --------------------------------------------------

    return {
        "occupancy_fraction": occupancy,
        "fields": mapped_fields,
        "x_edges": x_edges,
        "y_edges": y_edges,
    }