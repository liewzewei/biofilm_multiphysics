import numpy as np


def evaluate_mesh_quality(points, triangles):
    """
    Calculate quality metrics for every triangular element.

    points:    (n_nodes, 2)
    triangles: (n_triangles, 3)
    """

    # Coordinates of triangle vertices
    P = points[triangles]

    A = P[:, 0]
    B = P[:, 1]
    C = P[:, 2]

    # Edge vectors
    AB = B - A
    BC = C - B
    CA = A - C

    # Edge lengths
    a = np.linalg.norm(BC, axis=1)
    b = np.linalg.norm(CA, axis=1)
    c = np.linalg.norm(AB, axis=1)

    # Signed triangle area
    signed_area = 0.5 * (
        AB[:, 0] * (C[:, 1] - A[:, 1])
        - AB[:, 1] * (C[:, 0] - A[:, 0])
    )

    # Calculate angles using the cosine rule
    def angle(opposite, side1, side2):
        cos_theta = (
            side1**2 + side2**2 - opposite**2
        ) / (2 * side1 * side2)

        return np.degrees(
            np.arccos(np.clip(cos_theta, -1, 1))
        )

    angle_A = angle(a, b, c)
    angle_B = angle(b, a, c)
    angle_C = angle(c, a, b)

    min_angles = np.minimum.reduce([
        angle_A, angle_B, angle_C
    ])

    # Longest edge divided by shortest edge
    edges = np.stack([a, b, c], axis=1)

    aspect_ratio = (
        edges.max(axis=1) / edges.min(axis=1)
    )

    return {
        "signed_area": signed_area,
        "min_angle": min_angles,
        "aspect_ratio": aspect_ratio
    }