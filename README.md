# Vortex Panel Method with Boundary Layer Coupling

2D inviscid-viscous flow solver for airfoil aerodynamics, built as a
portfolio project.

## Physics
- Constant-strength vortex panel method (Katz & Plotkin §3.14)
- Thwaites laminar boundary layer (Thwaites 1949)
- Michel transition criterion (Michel 1951)
- Head turbulent boundary layer (Head 1958)
- Weak viscous-inviscid coupling via displacement thickness

## Validation
NACA 0012 against Abbott & von Doenhoff (1959) experimental data and XFOIL.

## Status
- [ ] Phase 1: Vortex panel method
- [ ] Phase 2: Boundary layer coupling
- [ ] Phase 3: Package cleanup and documentation