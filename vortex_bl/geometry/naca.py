# vortex_bl/geometry/naca.py
"""
NACA 4-digit airfoil geometry generator.

Theory: Abbott & von Doenhoff, "Theory of Wing Sections," Dover (1959), §6.
Thickness form: eq. 6.1, p. 113.
Panel spacing: cosine clustering per Katz & Plotkin, "Low Speed Aerodynamics,"
Cambridge (2001), §3.14.
"""
import numpy as np


def naca4_thickness(x: np.ndarray, thickness: float) -> np.ndarray:
    t = thickness
    yt = (t / 0.2) * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )
    return yt


def naca4_camber(x: np.ndarray, max_camber: float, max_camber_pos: float):
    m, p = max_camber, max_camber_pos
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    if m == 0.0 or p == 0.0:
        return yc, dyc_dx

    fwd = x < p
    aft = ~fwd

    yc[fwd] = (m / p**2) * (2 * p * x[fwd] - x[fwd]**2)
    yc[aft] = (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * x[aft] - x[aft]**2)

    dyc_dx[fwd] = (2 * m / p**2) * (p - x[fwd])
    dyc_dx[aft] = (2 * m / (1 - p)**2) * (p - x[aft])

    return yc, dyc_dx


def naca4_points(ndigits: str, n_panels: int = 100) -> tuple[np.ndarray, np.ndarray]:
    if len(ndigits) != 4 or not ndigits.isdigit():
        raise ValueError(f"Expected 4-digit NACA string, got '{ndigits}'")
    if n_panels % 2 != 0:
        raise ValueError("n_panels must be even.")

    max_camber = int(ndigits[0]) / 100.0
    max_camber_pos = int(ndigits[1]) / 10.0
    thickness = int(ndigits[2:]) / 100.0

    n_half = n_panels // 2

    theta = np.linspace(0, np.pi, n_half + 1)
    x_upper = 0.5 * (1 - np.cos(theta))
    x_upper = x_upper[::-1]  # reverse: TE -> LE

    yt_upper = naca4_thickness(x_upper, thickness)
    yc_upper, dyc_dx_upper = naca4_camber(x_upper, max_camber, max_camber_pos)

    theta_c = np.arctan(dyc_dx_upper)
    xu = x_upper - yt_upper * np.sin(theta_c)
    yu = yc_upper + yt_upper * np.cos(theta_c)

    x_lower = x_upper[::-1]
    yt_lower = naca4_thickness(x_lower, thickness)
    yc_lower, dyc_dx_lower = naca4_camber(x_lower, max_camber, max_camber_pos)
    theta_c_lower = np.arctan(dyc_dx_lower)
    xl = x_lower + yt_lower * np.sin(theta_c_lower)
    yl = yc_lower - yt_lower * np.cos(theta_c_lower)

    x = np.concatenate([xu, xl[1:]])
    y = np.concatenate([yu, yl[1:]])
    y[-1] = 0
    y[0] = 0

    return x, y