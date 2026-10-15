import numpy as np


def random_height_profile(
    x,
    mean_height=80.0,
    roughness=15.0,
    n_features=10,
    min_height=5.0,
    rng=None
):
    """
    Generate a continuous stochastic biofilm height profile.

    All geometric quantities use the same length units as x.
    """

    if rng is None:
        rng = np.random.default_rng()

    for attempt in range(1000):

        r = np.zeros_like(x, dtype=float)

        for _ in range(n_features):
            amplitude = rng.uniform(-1, 1)
            centre = rng.uniform(x.min(), x.max())
            spread = rng.uniform(15, 60)

            r += amplitude * np.exp(
                -0.5 * ((x - centre) / spread) ** 2
            )

        # Remove mean offset
        r -= np.trapezoid(r, x) / (x[-1] - x[0])

        # Normalize RMS roughness
        current_roughness = np.sqrt(
            np.trapezoid(r**2, x) / (x[-1] - x[0])
        )

        if current_roughness < 1e-12:
            continue

        r *= roughness / current_roughness

        h = mean_height + r

        if np.all(h >= min_height):
            return h

    raise RuntimeError("Could not generate a valid profile.")