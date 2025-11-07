import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Q4 shape functions
# ------------------------------------------------------------
def shape_Q4(xi, eta):
    # Shape functions
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta)
    ])

    # Derivatives wrt xi
    dN_dxi = 0.25 * np.array([
        -(1 - eta),
         (1 - eta),
         (1 + eta),
        -(1 + eta)
    ])

    # Derivatives wrt eta
    dN_deta = 0.25 * np.array([
        -(1 - xi),
        -(1 + xi),
         (1 + xi),
         (1 - xi)
    ])

    return N, dN_dxi, dN_deta


# ------------------------------------------------------------
# MAIN FEM ROUTINE
# ------------------------------------------------------------
def hw7_single_Q4_shear():
    # Inputs
    E = 200e9
    nu = 0.30
    L = 5.0
    h = 1.0
    b = 1.0
    tau = 200e3  # shear traction σ_xy

    # Plane stress constitutive matrix
    D = (E / (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0,  0, (1 - nu) / 2]
    ])

    # Node coordinates
    X = np.array([
        [0, 0],
        [L, 0],
        [L, h],
        [0, h]
    ], dtype=float)

    ndof = 8
    K = np.zeros((ndof, ndof))
    f = np.zeros(ndof)

    # 2×2 Gauss points
    gp = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    w  = [1, 1]

    # --------------------------------------------------------
    # ELEMENT STIFFNESS
    # --------------------------------------------------------
    for i in range(2):
        xi = gp[i]
        wi = w[i]
        for j in range(2):
            eta = gp[j]
            wj = w[j]

            N, dN_dxi, dN_deta = shape_Q4(xi, eta)

            # Jacobian
            dNparent = np.vstack((dN_dxi, dN_deta)).T  # 4×2
            J = dNparent.T @ X                         # 2×2
            detJ = np.linalg.det(J)
            if detJ <= 0:
                raise ValueError("Non-positive detJ")

            invJ = np.linalg.inv(J)

            # Gradients wrt x,y
            dN = dNparent @ invJ.T   # 4×2

            # Build B matrix
            B = np.zeros((3, 8))
            for a in range(4):
                B[0, 2*a]   = dN[a, 0]
                B[1, 2*a+1] = dN[a, 1]
                B[2, 2*a]   = dN[a, 1]
                B[2, 2*a+1] = dN[a, 0]

            # Stiffness accumulation
            K += (B.T @ D @ B) * b * detJ * wi * wj

    # --------------------------------------------------------
    # RIGHT EDGE TRACTION
    # --------------------------------------------------------
    tvec = np.array([0, tau])  # traction vector

    for j in range(2):
        eta = gp[j]
        wj = w[j]
        xi = 1.0

        N, _, dN_deta = shape_Q4(xi, eta)

        # Edge metric: derivative wrt eta → physical space
        x_eta = np.array([
            np.sum(dN_deta * X[:, 0]),
            np.sum(dN_deta * X[:, 1])
        ])
        ds_deta = np.linalg.norm(x_eta)

        # Consistent loads
        for a in range(4):
            idx = slice(2*a, 2*a+2)
            f[idx] += N[a] * tvec * b * ds_deta * wj

    # --------------------------------------------------------
    # APPLY BOUNDARY CONDITIONS
    # --------------------------------------------------------
    fixed = [0, 1, 6, 7]  # MATLAB 1-based -> Python 0-based
    free = [i for i in range(ndof) if i not in fixed]

    u = np.zeros(ndof)
    Kff = K[np.ix_(free, free)]
    ff = f[free]

    u[free] = np.linalg.solve(Kff, ff)

    # Reshape nodal displacement matrix
    U = u.reshape(4, 2)

    print("\nNodal Displacements (meters):")
    for i in range(4):
        print(f"Node {i+1}: ux = {U[i,0]:+.6e},   uy = {U[i,1]:+.6e}")

    tip_y = 0.5*(U[1,1] + U[2,1])
    print(f"\nTip deflection (avg of nodes 2 & 3): {tip_y:.6e} m")

    # --------------------------------------------------------
    # PLOT DEFORMED SHAPE
    # --------------------------------------------------------
    scale = 2000
    XY_def = X + scale * U

    plt.figure()
    plt.plot([X[:,0], X[:,0][[1,2,3,0]]], [X[:,1], X[:,1][[1,2,3,0]]],
             'k-', label="Original")
    plt.plot([XY_def[:,0], XY_def[:,0][[1,2,3,0]]], 
             [XY_def[:,1], XY_def[:,1][[1,2,3,0]]],
             'r--', label="Deformed (scaled)")
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    plt.title(f"Q4 Beam Under Shear (tau = {tau/1000:.0f} kPa)")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()


# Run
if __name__ == "__main__":
    hw7_single_Q4_shear()
