import numpy as np
from .panels import PanelGeometry


def solve_vortex_panels(panels: PanelGeometry, alpha_deg: float) -> dict:
    """
    Linear-strength vortex panel method.

    Direct translation of Katz & Plotkin, "Low Speed Aerodynamics,"
    2nd ed., Appendix D, Program 7 (linear-strength vortex, Neumann BC).

    M panels, N=M+1 nodes/unknowns.
    A is (N x N) influence matrix (normal BC + Kutta condition).
    B is (M x N) tangential influence matrix for velocity/Cp.
    """
    alpha = np.deg2rad(alpha_deg)

    # Program 7 reverses input to CW — we do the same
    XB = panels.x[::-1].copy()
    YB = panels.y[::-1].copy()

    M = panels.n_panels
    N = M + 1

    # Panel endpoints
    PT1x = XB[:-1];  PT1z = YB[:-1]
    PT2x = XB[1:];   PT2z = YB[1:]

    # Panel angles, control points, lengths
    DX = PT2x - PT1x;  DZ = PT2z - PT1z
    TH = np.arctan2(DZ, DX)
    CO1 = 0.5*(PT1x + PT2x)
    CO2 = 0.5*(PT1z + PT2z)
    DL  = np.sqrt(DX**2 + DZ**2)

    # A: N x N (M flow tangency rows + 1 Kutta row)
    # b: N RHS vector
    # B: M x N tangential influence
    A = np.zeros((N, N))
    b = np.zeros(N)
    B = np.zeros((M, N))

    for i in range(M):
        HOLDA = 0.0
        HOLDB = 0.0
        for j in range(M):
            XT  = CO1[i] - PT1x[j]
            ZT  = CO2[i] - PT1z[j]
            X2T = PT2x[j] - PT1x[j]
            Z2T = PT2z[j] - PT1z[j]

            X  =  XT*np.cos(TH[j]) + ZT*np.sin(TH[j])
            Z  = -XT*np.sin(TH[j]) + ZT*np.cos(TH[j])
            X2 =  X2T*np.cos(TH[j]) + Z2T*np.sin(TH[j])

            # Save panel lengths (on first pass i=0)
            if i == 0:
                DL[j] = X2

            R1  = np.sqrt(X**2 + Z**2)
            R2  = np.sqrt((X-X2)**2 + Z**2)
            TH1 = np.arctan2(Z, X)
            TH2 = np.arctan2(Z, X-X2)

            if i == j:
                U1L = -0.5*(X-X2)/X2
                U2L =  0.5*X/X2
                W1L = -0.15916
                W2L =  0.15916
            else:
                U1L = -(Z*np.log(R2/R1) + X*(TH2-TH1) - X2*(TH2-TH1)) \
                       / (6.28319*X2)
                U2L =  (Z*np.log(R2/R1) + X*(TH2-TH1)) \
                       / (6.28319*X2)
                W1L = -((X2 - Z*(TH2-TH1)) - X*np.log(R1/R2)
                       + X2*np.log(R1/R2)) / (6.28319*X2)
                W2L =  ((X2 - Z*(TH2-TH1)) - X*np.log(R1/R2)) \
                       / (6.28319*X2)

            # Rotate to global frame
            U1 =  U1L*np.cos(-TH[j]) + W1L*np.sin(-TH[j])
            U2 =  U2L*np.cos(-TH[j]) + W2L*np.sin(-TH[j])
            W1 = -U1L*np.sin(-TH[j]) + W1L*np.cos(-TH[j])
            W2 = -U2L*np.sin(-TH[j]) + W2L*np.cos(-TH[j])

            # Assemble A and B
            if j == 0:
                A[i, 0]  = -U1*np.sin(TH[i]) + W1*np.cos(TH[i])
                HOLDA    = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
                B[i, 0]  =  U1*np.cos(TH[i]) + W1*np.sin(TH[i])
                HOLDB    =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])
            elif j == M-1:
                A[i, M-1] = -U1*np.sin(TH[i]) + W1*np.cos(TH[i]) + HOLDA
                A[i, M]   = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
                B[i, M-1] =  U1*np.cos(TH[i]) + W1*np.sin(TH[i]) + HOLDB
                B[i, M]   =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])
            else:
                A[i, j]  = -U1*np.sin(TH[i]) + W1*np.cos(TH[i]) + HOLDA
                HOLDA    = -U2*np.sin(TH[i]) + W2*np.cos(TH[i])
                B[i, j]  =  U1*np.cos(TH[i]) + W1*np.sin(TH[i]) + HOLDB
                HOLDB    =  U2*np.cos(TH[i]) + W2*np.sin(TH[i])

        # RHS for row i
        b[i] = np.cos(alpha)*np.sin(TH[i]) - np.sin(alpha)*np.cos(TH[i])

    # Kutta condition: G[0] + G[N-1] = 0
    A[M, 0]   = 1.0
    A[M, N-1] = 1.0
    b[M]      = 0.0

    # Solve
    G = np.linalg.solve(A, b)

    # Surface velocity and Cp
    VEL = np.zeros(M)
    Cp  = np.zeros(M)
    CL  = 0.0
    for i in range(M):
        vel = 0.0
        for j in range(N):
            vel += B[i, j] * G[j]
        vel += np.cos(alpha)*np.cos(TH[i]) + np.sin(alpha)*np.sin(TH[i])
        VEL[i] = vel
        Cp[i]  = 1.0 - vel**2
        CL    += vel * DL[i]

    # Reverse output back to original CCW ordering
    VEL = VEL[::-1]
    Cp  = Cp[::-1]

    return {
        'gamma': G,
        'Cp':    Cp,
        'Vt':    VEL,
        'Cl':    2.0 * CL,
        'alpha': alpha_deg,
    }