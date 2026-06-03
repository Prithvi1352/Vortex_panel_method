# vortex_bl/postprocess/plots.py
import numpy as np
import matplotlib.pyplot as plt
from ..panel.panels import PanelGeometry


def plot_cp(panels: PanelGeometry, result: dict, ax=None, label=None, **kwargs):
    """
    Plot Cp distribution against x/c.
    Inverted y-axis (suction up) per aerodynamic convention.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    ax.plot(panels.xc, result['Cp'], label=label, **kwargs)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel('x/c')
    ax.set_ylabel(r'$C_p$')
    ax.set_title(f"Cp distribution — "
                 f"α = {result['alpha']:.1f}°, "
                 f"$C_l$ = {result['Cl']:.4f}")
    ax.axhline(1.0, color='k', linewidth=0.5, linestyle='--', alpha=0.4)
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_polar(alphas, cls, ax=None, label=None, **kwargs):
    """
    Plot Cl vs alpha polar.
    Overlays thin airfoil theory reference.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    ax.plot(alphas, cls, label=label, **kwargs)

    # Thin airfoil theory reference
    alpha_ref = np.linspace(alphas[0], alphas[-1], 200)
    cl_thin = 2 * np.pi * np.sin(np.deg2rad(alpha_ref))
    ax.plot(alpha_ref, cl_thin, 'k--', alpha=0.5,
            label=r'Thin airfoil: $2\pi\sin\alpha$')

    ax.set_xlabel(r'$\alpha$ (deg)')
    ax.set_ylabel(r'$C_l$')
    ax.set_title('Lift polar — NACA 0012 inviscid')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


def plot_airfoil(panels: PanelGeometry, ax=None):
    """Plot airfoil geometry."""
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 3))

    ax.plot(panels.x, panels.y, 'k-', linewidth=1.5)
    ax.plot(panels.xc, panels.yc, 'r.', markersize=3, label='Control points')
    ax.set_aspect('equal')
    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.set_title('NACA 0012 — panel geometry')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    return ax