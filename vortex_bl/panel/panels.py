# vortex_bl/panel/panels.py
"""
Panel geometry from airfoil node coordinates.

Follows the panel conventions in Katz & Plotkin, "Low Speed Aerodynamics,"
Cambridge (2001), §3.14, and Moran, "Introduction to Theoretical and
Computational Aerodynamics," Wiley (1984), §5.3.
"""
import numpy as np
from dataclasses import dataclass


@dataclass
class PanelGeometry:
    x: np.ndarray
    y: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    length: np.ndarray
    cos_phi: np.ndarray
    sin_phi: np.ndarray
    nx: np.ndarray
    ny: np.ndarray
    phi: np.ndarray
    n_panels: int


def build_panels(x_nodes: np.ndarray, y_nodes: np.ndarray) -> PanelGeometry:
    n = len(x_nodes) - 1

    dx = np.diff(x_nodes)
    dy = np.diff(y_nodes)

    length = np.hypot(dx, dy)
    cos_phi = dx / length
    sin_phi = dy / length

    xc = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    yc = 0.5 * (y_nodes[:-1] + y_nodes[1:])

    # Outward normal: rotate tangent 90° counter-clockwise
    # (Katz & Plotkin eq. 3.76)
    nx = sin_phi
    ny = -cos_phi

    phi = np.arctan2(sin_phi, cos_phi)

    return PanelGeometry(
        x=x_nodes, y=y_nodes,
        xc=xc, yc=yc,
        length=length,
        cos_phi=cos_phi, sin_phi=sin_phi,
        nx=nx, ny=ny,
        phi=phi,
        n_panels=n,
    )