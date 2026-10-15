"""
Estimate biofilm-substrate interfacial tractions.

Uses element-averaged stresses from triangular FEM elements.

Units:
    Coordinates: micrometres
    Stress: Pa
    Traction: Pa
"""

import numpy as np


def calculate_interface_traction(
    points,
    triangles,
    stress,
    substrate_nodes,
):
    """
    Calculate normal and tangential traction along the
    biofilm-substrate interface.

    Parameters
    ----------
    points : ndarray (N, 2)
        FEM mesh coordinates in micrometres.

    triangles : ndarray (M, 3)
        Triangle connectivity.

    stress : ndarray (M, 3)
        Element-averaged stress components:
        [sigma_xx, sigma_yy, tau_xy].

    substrate_nodes : ndarray
        Ordered substrate boundary node indices.

    Returns
    -------
    dict
        Segment midpoint coordinates, lengths,
        normal traction, and shear traction.
    """

    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int32)
    stress = np.asarray(stress, dtype=float)
    substrate_nodes = np.asarray(substrate_nodes, dtype=np.int32)

    if stress.shape != (len(triangles), 3):
        raise ValueError("Stress must have shape (N_elements, 3).")

    # --------------------------------------------------
    # 1. Construct edge-to-element lookup
    # --------------------------------------------------

    edge_to_element = {}

    for element_id, tri in enumerate(triangles):

        edges = [
            (tri[0], tri[1]),
            (tri[1], tri[2]),
            (tri[2], tri[0]),
        ]

        for a, b in edges:

            edge = tuple(sorted((int(a), int(b))))

            edge_to_element.setdefault(
                edge, []
            ).append(element_id)

    # --------------------------------------------------
    # 2. Calculate interface tractions
    # --------------------------------------------------

    midpoints = []
    lengths = []
    normal_tractions = []
    shear_tractions = []

    for i in range(len(substrate_nodes) - 1):

        a = int(substrate_nodes[i])
        b = int(substrate_nodes[i + 1])

        edge = tuple(sorted((a, b)))

        adjacent_elements = edge_to_element.get(edge, [])

        if len(adjacent_elements) != 1:
            raise ValueError(
                f"Expected one adjacent element for boundary edge {edge}; "
                f"found {len(adjacent_elements)}."
            )

        element_id = adjacent_elements[0]

        # Boundary segment geometry
        pa = points[a]
        pb = points[b]

        edge_vector = pb - pa

        length = np.linalg.norm(edge_vector)

        if length <= 0:
            raise ValueError("Zero-length boundary edge.")

        # Unit tangent
        tangent = edge_vector / length

        # Unit normal pointing into the biofilm
        normal = np.array([
            -tangent[1],
            tangent[0],
        ])

        # Element stress tensor
        sxx, syy, txy = stress[element_id]

        sigma = np.array([
            [sxx, txy],
            [txy, syy],
        ])

        # Cauchy traction
        traction = sigma @ normal

        # Normal component
        tn = normal @ traction

        # Tangential component
        ts = tangent @ traction

        midpoints.append((pa + pb) / 2)
        lengths.append(length)
        normal_tractions.append(tn)
        shear_tractions.append(ts)

    return {
        "midpoints": np.asarray(midpoints),
        "lengths": np.asarray(lengths),
        "normal_traction": np.asarray(normal_tractions),
        "shear_traction": np.asarray(shear_tractions),
    }