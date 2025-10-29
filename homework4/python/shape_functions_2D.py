# ========================================================================
# ENGR 4350/6350 — 2D Isoparametric Shape Functions (Q4 / Q9)
# Trey Brown — Fall 2025
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# N2D — Shape Functions
# ----------------------------------------------------------------------
def N2D(xi: float, eta: float, nn: int) -> np.ndarray:
    """
    Return 2D shape functions for Q4 or Q9 elements.
    Output: 1×nn numpy array of shape function values.
    """
    if nn == 4:
        N = 0.25 * np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ])
    elif nn == 9:
        Lxi  = np.array([0.5 * xi * (xi - 1), 1 - xi**2, 0.5 * xi * (xi + 1)])
        Leta = np.array([0.5 * eta * (eta - 1), 1 - eta**2, 0.5 * eta * (eta + 1)])

        N = np.array([
            Lxi[0]*Leta[0], Lxi[2]*Leta[0], Lxi[2]*Leta[2], Lxi[0]*Leta[2],
            Lxi[1]*Leta[0], Lxi[2]*Leta[1], Lxi[1]*Leta[2], Lxi[0]*Leta[1],
            Lxi[1]*Leta[1]
        ])
    else:
        raise ValueError("N2D: nn must be 4 or 9")
    return N


# ----------------------------------------------------------------------
# GN2D — Derivatives of Shape Functions
# ----------------------------------------------------------------------
def GN2D(xi: float, eta: float, nn: int) -> np.ndarray:
    """
    Return derivatives of shape functions wrt (xi, eta).
    Output: 2×nn matrix [dN/dxi; dN/deta].
    """
    if nn == 4:
        dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
        GN = np.vstack((dN_dxi, dN_deta))

    elif nn == 9:
        # 1D bases and derivatives
        Lxi   = np.array([0.5*xi*(xi-1), 1 - xi**2, 0.5*xi*(xi+1)])
        dLxi  = np.array([xi - 0.5, -2*xi, xi + 0.5])
        Leta  = np.array([0.5*eta*(eta-1), 1 - eta**2, 0.5*eta*(eta+1)])
        dLeta = np.array([eta - 0.5, -2*eta, eta + 0.5])

        dN_dxi = np.array([
            dLxi[0]*Leta[0], dLxi[2]*Leta[0], dLxi[2]*Leta[2], dLxi[0]*Leta[2],
            dLxi[1]*Leta[0], dLxi[2]*Leta[1], dLxi[1]*Leta[2], dLxi[0]*Leta[1],
            dLxi[1]*Leta[1]
        ])

        dN_deta = np.array([
            Lxi[0]*dLeta[0], Lxi[2]*dLeta[0], Lxi[2]*dLeta[2], Lxi[0]*dLeta[2],
            Lxi[1]*dLeta[0], Lxi[2]*dLeta[1], Lxi[1]*dLeta[2], Lxi[0]*dLeta[1],
            Lxi[1]*dLeta[1]
        ])

        GN = np.vstack((dN_dxi, dN_deta))

    else:
        raise ValueError("GN2D: nn must be 4 or 9")

    return GN


# ----------------------------------------------------------------------
# Plot Comparison — Visualize 4- and 9-node shape functions
# ----------------------------------------------------------------------
def plot_shape_functions_compare():
    xiplot = np.linspace(-1, 1, 25)
    etaplot = np.linspace(-1, 1, 25)
    Xi, Eta = np.meshgrid(xiplot, etaplot)

    for nn in [4, 9]:
        theta = np.zeros_like(Xi)
        dth_dxi = np.zeros_like(Xi)
        dth_deta = np.zeros_like(Xi)
        d = np.arange(1, nn + 1)  # arbitrary nodal values (for visualization)

        for i in range(Xi.size):
            xi = Xi.flat[i]
            eta = Eta.flat[i]
            N = N2D(xi, eta, nn)
            GN = GN2D(xi, eta, nn)
            theta.flat[i] = N @ d
            dth_dxi.flat[i] = GN[0, :] @ d
            dth_deta.flat[i] = GN[1, :] @ d

        # --- Surface plot ---
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(Xi, Eta, theta, cmap='viridis', edgecolor='none')
        ax.set_title(f'Shape Function Surface ({nn}-node)')
        ax.set_xlabel('ξ')
        ax.set_ylabel('η')
        ax.set_zlabel('θ')
        plt.tight_layout()
        plt.show()

        # --- Gradient field plot ---
        plt.figure()
        plt.quiver(Xi, Eta, dth_dxi, dth_deta)
        plt.title(f'Gradient Field ({nn}-node)')
        plt.xlabel('ξ')
        plt.ylabel('η')
        plt.axis('equal')
        plt.grid(True)
        plt.show()


# ----------------------------------------------------------------------
# Run comparison if executed directly
# ----------------------------------------------------------------------
if __name__ == "__main__":
    plot_shape_functions_compare()
