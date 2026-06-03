# run_day2.py
"""
Day 2: Cp distribution, lift polar, geometry plot.
Run from project root: python run_day2.py
"""
import numpy as np
import matplotlib.pyplot as plt

from vortex_bl.geometry.naca import naca4_points
from vortex_bl.panel.panels import build_panels
from vortex_bl.panel.solver import solve_vortex_panels
from vortex_bl.postprocess.plots import plot_cp, plot_polar, plot_airfoil


def compute_polar(alpha_range, n_panels=100):
    x, y = naca4_points('0012', n_panels=n_panels)
    panels = build_panels(x, y)
    alphas, cls = [], []
    for alpha in alpha_range:
        r = solve_vortex_panels(panels, alpha)
        alphas.append(alpha)
        cls.append(r['Cl'])
        print(f"  alpha={alpha:+6.1f}  Cl={r['Cl']:+.4f}")
    return np.array(alphas), np.array(cls)


if __name__ == '__main__':

    # ── 1. Airfoil geometry ───────────────────────────────────────────
    x, y = naca4_points('0012', n_panels=100)
    panels = build_panels(x, y)
    fig, ax = plt.subplots(figsize=(10, 3))
    plot_airfoil(panels, ax=ax)
    plt.tight_layout()
    plt.savefig('output/airfoil_geometry.png', dpi=150)
    plt.close()
    print('Saved: output/airfoil_geometry.png')

    # ── 2. Cp at multiple angles ──────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, alpha in zip(axes, [0, 5, 10]):
        r = solve_vortex_panels(panels, float(alpha))
        plot_cp(panels, r, ax=ax, color='steelblue')
    plt.tight_layout()
    plt.savefig('output/cp_distributions.png', dpi=150)
    plt.close()
    print('Saved: output/cp_distributions.png')

    # ── 3. Cp convergence with N ──────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['lightblue', 'steelblue', 'navy', 'darkred']
    for n, color in zip([20, 50, 100, 200], colors):
        x, y = naca4_points('0012', n_panels=n)
        p = build_panels(x, y)
        r = solve_vortex_panels(p, 5.0)
        plot_cp(p, r, ax=ax, color=color, label=f'N={n}')
    ax.set_title('Cp convergence — NACA 0012, α=5°')
    plt.tight_layout()
    plt.savefig('output/cp_convergence.png', dpi=150)
    plt.close()
    print('Saved: output/cp_convergence.png')

    # ── 4. Cl polar ───────────────────────────────────────────────────
    print('\nComputing polar...')
    alpha_range = np.arange(-10, 16, 1)
    alphas, cls = compute_polar(alpha_range, n_panels=100)

    fig, ax = plt.subplots(figsize=(6, 6))
    plot_polar(alphas, cls, ax=ax, color='steelblue',
               marker='o', markersize=4, label='VPM N=100')
    plt.tight_layout()
    plt.savefig('output/cl_polar.png', dpi=150)
    plt.close()
    print('Saved: output/cl_polar.png')

    # ── 5. Summary ────────────────────────────────────────────────────
    x, y = naca4_points('0012', n_panels=100)
    panels = build_panels(x, y)
    r5 = solve_vortex_panels(panels, 5.0)
    dcl = r5['Cl'] / np.deg2rad(5.0)
    print(f'\nLift slope dCl/dalpha = {dcl:.4f} rad⁻¹')
    print(f'Thin airfoil theory   = {2*np.pi:.4f} rad⁻¹')
    print(f'Error                 = {abs(dcl - 2*np.pi)/( 2*np.pi)*100:.1f}%')