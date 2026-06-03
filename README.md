# 2D Vortex Panel Method with Boundary Layer Coupling

A from-scratch implementation of a 2D inviscid-viscous aerodynamic solver for airfoil analysis, built as an independent portfolio project to demonstrate physics understanding and scientific computing ability.

---

## What This Project Does

This solver computes the pressure distribution, lift coefficient, and drag coefficient of a 2D airfoil in low-speed flow. It combines an inviscid vortex panel method with a viscous boundary layer solver, producing results that account for skin friction drag and viscous displacement effects.

**Outputs:**
- Cp distribution over the airfoil surface
- Cl vs angle of attack polar (inviscid and viscous)
- Cd from boundary layer integral methods
- Boundary layer thickness distribution δ, δ*, θ, H along the chord
- Transition location prediction

**Validated against** XFOIL (inviscid mode) and Abbott & von Doenhoff (1959) experimental data for NACA 0012.

---

## Physics and Methods

### Phase 1 — Inviscid Vortex Panel Method

The inviscid solver is based on the **linear-strength vortex panel method** from Katz & Plotkin, *Low Speed Aerodynamics* (2nd ed., Cambridge 2001), Appendix D, Program 7.

**Governing equation:** Laplace's equation for the velocity potential ∇²φ = 0, solved by distributing vortex singularities on the airfoil surface.

**Discretisation:** The airfoil is divided into M panels. The vortex sheet strength γ varies linearly between node values, giving N = M+1 unknowns. The no-penetration boundary condition (zero normal velocity at each panel midpoint) provides M equations. The Kutta condition γ₁ + γ_N = 0 provides the (M+1)th equation, enforcing smooth flow at the trailing edge.

**Influence coefficients:** For each panel j, the velocity induced at control point i is computed analytically by integrating the 2D Biot-Savart law along the panel. The linear vortex distribution splits into start-node and end-node contributions (U1L, W1L) and (U2L, W2L), rotated to global frame and assembled into the influence matrices A (normal) and B (tangential).

**Surface velocity and pressure:** After solving for Γ, the tangential velocity at each panel midpoint is:

```
Vt_i = Σ B_ij * Γ_j + cos(α)cos(φᵢ) + sin(α)sin(φᵢ)
Cp_i = 1 - Vt_i²
```

**Lift coefficient:** Cl = 2 Σ Vt_i Δsᵢ (Kutta-Joukowski, unit chord, unit freestream)

**Why linear-strength and not constant-strength:** The constant-strength vortex panel has a zero diagonal in the influence matrix and a rank deficiency when the Kutta condition is applied. The linear-strength method has a nonzero self-induced tangential velocity (±0.15916) giving a well-conditioned diagonal regardless of panel spacing.

### Phase 2 — Boundary Layer Coupling (in progress)

The boundary layer solver uses the von Kármán integral momentum equation as the foundation, with separate methods for the laminar and turbulent regions.

**Laminar boundary layer — Thwaites (1949):**

The Thwaites method solves the von Kármán equation via a universal correlation. The momentum thickness evolves as:

```
θ²(x) = (0.45ν / Ue⁶) ∫₀ˣ Ue⁵ dx
```

The shape factor H and skin friction Cf are recovered from tabulated correlations in λ = (θ²/ν)(dUe/dx).

**Transition — Michel (1951):**

Transition from laminar to turbulent is predicted when:

```
Re_θ ≥ 2.9 Re_x^0.4
```

where Re_θ = Ue θ/ν and Re_x = Ue x/ν.


**Turbulent boundary layer — Head (1958):**

Head's entrainment method adds a second ODE for the shape factor H₁, governing the rate at which irrotational fluid is entrained into the turbulent region. The coupled system is integrated from the transition point to the trailing edge.


**Viscous-inviscid coupling:** Weak coupling via displacement thickness. The inviscid surface velocity distribution is modified by the displacement effect δ* of the boundary layer, and the panel method is re-run iteratively until convergence.

---


## Usage

```python
from vortex_bl.geometry.naca import naca4_points
from vortex_bl.panel.panels import build_panels
from vortex_bl.panel.solver import solve_vortex_panels

# Generate NACA 0012 geometry
x, y = naca4_points('0012', n_panels=100)
panels = build_panels(x, y)

# Solve at alpha = 5 degrees
result = solve_vortex_panels(panels, alpha_deg=5.0)

print(f"Cl = {result['Cl']:.4f}")
# Cl = 0.5999

# Access full distributions
xc = panels.xc       # control point x/c locations
Cp = result['Cp']    # pressure coefficient
Vt = result['Vt']    # surface tangential velocity
```

## References

**Panel method theory:**
- Katz, J. & Plotkin, A. (2001). *Low Speed Aerodynamics* (2nd ed.). Cambridge University Press. — Primary reference for the linear-strength vortex panel formulation (§11.2, Program 7).
- Kuethe, A.M. & Chow, C.Y. (1998). *Foundations of Aerodynamics* (5th ed.). Wiley. — Kutta-Joukowski theorem, thin airfoil theory context.
- Moran, J. (1984). *Introduction to Theoretical and Computational Aerodynamics*. Wiley. — Panel method convergence and accuracy discussion.
- Hess, J.L. & Smith, A.M.O. (1967). Calculation of potential flow about arbitrary bodies. *Progress in Aerospace Sciences*, 8, 1–138. — Foundational panel method paper.

**Boundary layer:**
- Thwaites, B. (1949). Approximate calculation of the laminar boundary layer. *Aeronautical Quarterly*, 1, 245–280.
- Michel, R. (1951). *Étude de la Transition sur les Profils d'Aile*. ONERA Report 1/1578A.
- Head, M.R. (1958). *Entrainment in the Turbulent Boundary Layer*. ARC R&M 3152.
- White, F.M. (2005). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill. — Modern treatment of Thwaites and Head methods.
- Drela, M. (2014). *Flight Vehicle Aerodynamics*. MIT Press. — Viscous-inviscid coupling methodology.

**Validation data:**
- Abbott, I.H. & von Doenhoff, A.E. (1959). *Theory of Wing Sections*. Dover. — NACA 0012 experimental Cp and polar data (Appendix IV).
- Drela, M. XFOIL 6.99. web.mit.edu/drela/Public/web/xfoil — Inviscid and viscous reference results.

---

## Status

- [x] Phase 1: Linear-strength vortex panel method — complete
- [x] Validation against XFOIL inviscid — 0.57% agreement
- [ ] Phase 2: Boundary layer coupling — in progress
- [ ] Phase 3: Package cleanup, documentation, GitHub publication
