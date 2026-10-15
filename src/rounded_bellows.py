import numpy as np


def rounded_bellows_profile(
    segment_length,
    angle_deg,
    radius,
    n_straight=100,
    n_arc=101,
):
    """
    Generate a symmetric V-shaped bellows with a rounded apex.

    Parameters
    ----------
    segment_length : float
        Length of each straight segment (μm).

    angle_deg : float
        Opening angle of the V (degrees).

    radius : float
        Radius of curvature at the apex (μm).

    n_straight : int
        Number of sampling points per straight segment.

    n_arc : int
        Number of sampling points along the rounded apex.

    Returns
    -------
    points : ndarray, shape (N, 2)
        Bellows surface coordinates (μm).

    s : ndarray
        Cumulative arc-length coordinate (μm).
    """

    theta = np.deg2rad(angle_deg)

    if not (0 < theta < np.pi):
        raise ValueError("Angle must be between 0 and 180 degrees.")

    if segment_length <= 0 or radius <= 0:
        raise ValueError("Lengths and radius must be positive.")

    if n_straight < 2 or n_arc < 3:
        raise ValueError("Insufficient sampling points.")

    # Angle between each straight segment and horizontal
    phi = 0.5 * (np.pi - theta)

    # Circular arc tangent to both straight segments.
    # Arc turning angle equals pi - theta.
    arc_half_angle = phi

    # Place the circle centre at the origin.
    # The upper arc is parameterized from left to right.
    arc_angles = np.linspace(
        np.pi / 2 + arc_half_angle,
        np.pi / 2 - arc_half_angle,
        n_arc,
    )

    arc_x = radius * np.cos(arc_angles)
    arc_y = radius * np.sin(arc_angles)

    # Left tangent point
    left_tangent = np.array([arc_x[0], arc_y[0]])

    # Right tangent point
    right_tangent = np.array([arc_x[-1], arc_y[-1]])

    # Straight segment directions, from endpoints toward apex
    left_direction = np.array([
        np.cos(phi),
        np.sin(phi),
    ])

    right_direction = np.array([
        -np.cos(phi),
        np.sin(phi),
    ])

    # Construct left straight segment
    t_left = np.linspace(
        segment_length,
        0,
        n_straight,
        endpoint=False,
    )

    left_points = (
        left_tangent[None, :]
        - t_left[:, None] * left_direction[None, :]
    )

    # Construct right straight segment
    t_right = np.linspace(
        0,
        segment_length,
        n_straight,
    )

    right_points = (
        right_tangent[None, :]
        - t_right[:, None] * right_direction[None, :]
    )

    # Combine segments
    arc_points = np.column_stack([arc_x, arc_y])

    points = np.vstack([
        left_points,
        arc_points,
        right_points[1:],
    ])

    # Shift minimum y to zero
    points[:, 1] -= points[:, 1].min()

    # Calculate cumulative arc length
    ds = np.linalg.norm(
        np.diff(points, axis=0),
        axis=1,
    )

    s = np.concatenate([
        [0.0],
        np.cumsum(ds),
    ])

    return points, s

def deform_rounded_bellows(
    segment_length,
    initial_angle_deg,
    final_angle_deg,
    initial_radius,
    n_straight=100,
    n_arc=101,
):
    """
    Generate initial and final bellows configurations
    while preserving total material arc length.

    Assumptions:
        - Straight segment lengths remain constant.
        - Rounded fold remains circular.
        - Fold arc length remains constant.
        - Fold radius changes with opening angle.

    All lengths are in micrometres.
    """

    theta0 = np.deg2rad(initial_angle_deg)
    theta1 = np.deg2rad(final_angle_deg)

    if not (0 < theta0 < np.pi and 0 < theta1 < np.pi):
        raise ValueError(
            "Opening angles must be between 0 and 180 degrees."
        )

    # Preserve rounded fold arc length
    radius_final = (
        initial_radius
        * (np.pi - theta0)
        / (np.pi - theta1)
    )

    # Initial configuration
    points_initial, s_initial = rounded_bellows_profile(
        segment_length=segment_length,
        angle_deg=initial_angle_deg,
        radius=initial_radius,
        n_straight=n_straight,
        n_arc=n_arc,
    )

    # Final configuration
    points_final, s_final = rounded_bellows_profile(
        segment_length=segment_length,
        angle_deg=final_angle_deg,
        radius=radius_final,
        n_straight=n_straight,
        n_arc=n_arc,
    )

    return {
        "initial_points": points_initial,
        "final_points": points_final,
        "s_initial": s_initial,
        "s_final": s_final,
        "radius_initial": initial_radius,
        "radius_final": radius_final,
    }
    
def resample_bellows(points, s, s_new):
    """
    Resample a bellows profile onto specified material coordinates.

    Parameters
    ----------
    points : ndarray (N, 2)
        Bellows coordinates in micrometres.

    s : ndarray (N,)
        Existing cumulative material coordinates.

    s_new : ndarray (M,)
        Desired material-coordinate positions.

    Returns
    -------
    ndarray (M, 2)
        Resampled bellows coordinates.
    """

    points = np.asarray(points, dtype=float)
    s = np.asarray(s, dtype=float)
    s_new = np.asarray(s_new, dtype=float)

    x_new = np.interp(
        s_new,
        s,
        points[:, 0]
    )

    y_new = np.interp(
        s_new,
        s,
        points[:, 1]
    )

    return np.column_stack([x_new, y_new])