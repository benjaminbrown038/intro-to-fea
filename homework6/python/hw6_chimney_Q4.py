# ========================================================================
# ENGR 4350/6350 — Homework 6: Chimney Conduction (Q4 FEM)
# Python translation of hw6_chimney_Q4.m
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# ----------------------------
# Helper functions
# ----------------------------
def N_Q4(xi: float, eta: float):
    """Q4 shape functions and parent derivatives."""
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)
    N = np.array([N1, N2, N3, N4])

    dN_dxi = 0.25 * np.array([-(1 - eta), +(1 - eta), +(1 + eta), -(1 + eta)])  # (4,)
    dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), +(1 + xi), +(1 - xi)])    # (4,)
    return N, dN_dxi, dN_deta


def gauss2x2():
    g = 1.0 / np.sqrt(3.0)
    gp = np.array([[-g, -g],
                   [ g, -g],
                   [ g,  g],
                   [-g,  g]], dtype=float)   # (4,2) -> [xi, eta]
    gw = np.ones(4, dtype=float)            # (4,)
    return gp, gw


def trisurf_draped_as_contour(nodes: np.ndarray, conn: np.ndarray, u: np.ndarray):
    """
    Split each quad into two triangles and plot a filled triangulation with scalar u.
    nodes: (n,2), conn: (m,4)  [0-based], u: (n,)
    """
    # split each quad [1 2 3 4] -> triangles [1 2 3] and [1 3 4]
    tri = np.vstack([conn[:, [0, 1, 2]], conn[:, [0, 2, 3]]])
    triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tri)
    plt.tripcolor(triang, u, shading="gouraud", edgecolors='k', linewidth=0.4)
    plt.gca().set_aspect('equal', 'box')


# ----------------------------
# Main
# ----------------------------
def hw6_chimney_Q4():
    # Material properties
    k_brick = 0.9
    k_concrete = 2.0

    # Boundary temperatures (°C)
    T_in = 140.0   # inner arc (hot gas)
    T_out = 10.0   # outer arc (ambient)

    # Node coordinates (approx. 1/8 circular region)
    X = np.array([[0.30, 0.35, 0.40],
                  [0.30, 0.40, 0.50],
                  [0.30, 0.45, 0.60]], dtype=float)

    Y = np.array([[0.40, 0.40, 0.40],
                  [0.50, 0.50, 0.50],
                  [0.60, 0.60, 0.60]], dtype=float)

    nodes = np.column_stack([X.ravel(order='C'), Y.ravel(order='C')])  # (9,2)
    ndof = nodes.shape[0]

    # Element connectivity (counterclockwise) — convert to 0-based
    conn = np.array([
        [1, 4, 5, 2],
        [2, 5, 6, 3],
        [4, 7, 8, 5],
        [5, 8, 9, 6]
    ], dtype=int) - 1

    # Material assignment per element
    k_elem = np.array([k_brick, k_concrete, k_brick, k_concrete], dtype=float)

    # --- Boundary nodes on inner/outer arcs (diagonal ~ circular arcs) ---
    tol = 1e-12
    x = nodes[:, 0]
    y = nodes[:, 1]
    on_diag = np.abs(x - y) < tol

    inner_arc_nodes = np.where(on_diag & (x <= 0.40 + tol))[0]
    outer_arc_nodes = np.where(on_diag & (x >= 0.60 - tol))[0]

    BC_nodes = np.concatenate([inner_arc_nodes, outer_arc_nodes])
    BC_vals = np.concatenate([np.full(inner_arc_nodes.size, T_in),
                              np.full(outer_arc_nodes.size, T_out)])

    # --- Assembly ---
    K = np.zeros((ndof, ndof), dtype=float)
    f = np.zeros(ndof, dtype=float)

    gp, gw = gauss2x2()

    for e, kappa in enumerate(k_elem):
        en = conn[e, :]              # (4,)
        xe = nodes[en, 0]            # (4,)
        ye = nodes[en, 1]            # (4,)
        ke = np.zeros((4, 4), dtype=float)

        for ig in range(4):
            xi, eta = gp[ig, 0], gp[ig, 1]
            w = gw[ig]

            N, dN_dxi, dN_deta = N_Q4(xi, eta)

            # Jacobian J = [dN/dxi; dN/deta] * [x y]
            J = np.array([[np.dot(dN_dxi, xe), np.dot(dN_dxi, ye)],
                          [np.dot(dN_deta, xe), np.dot(dN_deta, ye)]], dtype=float)
            detJ = np.linalg.det(J)
            if detJ <= 0:
                raise ValueError(f"Non-positive detJ in element {e+1}")
            invJ = np.linalg.inv(J)

            # Derivatives wrt global coords: [dNdx dNdy] = [dN_dxi dN_deta] * invJ^T
            dN_parent = np.column_stack([dN_dxi, dN_deta])  # (4,2)
            dN_global = dN_parent @ invJ.T                 # (4,2)
            B = dN_global.T                                # (2,4)

            # ke += (B^T * kappa * I * B) * detJ * w
            ke += (B.T @ (kappa * np.eye(2)) @ B) * detJ * w

        # Scatter-add to global K
        K[np.ix_(en, en)] += ke

    # --- Apply Dirichlet BCs and solve ---
    free = np.ones(ndof, dtype=bool)
    free[BC_nodes] = False

    uc = np.full(ndof, np.nan, dtype=float)
    uc[BC_nodes] = BC_vals

    Kff = K[np.ix_(free, free)]
    Kfc = K[np.ix_(free, ~free)]
    rhs = f[free] - Kfc @ uc[~free]

    uf = np.linalg.solve(Kff, rhs)

    u = np.zeros(ndof, dtype=float)
    u[free] = uf
    u[~free] = uc[~free]

    # --- Output ---
    print("\nNodal Temperatures (°C):")
    for n in range(ndof):
        print(f"  Node {n+1:2d} ({nodes[n,0]:.3f}, {nodes[n,1]:.3f}):  T = {u[n]:8.4f}")

    print("\nNode numbering (x, y):")
    for n in range(ndof):
        print(f"  {n+1:2d}: ({nodes[n,0]:.3f}, {nodes[n,1]:.3f})")

    print("\nElement connectivity (Q4):")
    for e in range(conn.shape[0]):
        a, b, c, d = (conn[e, :] + 1)  # back to 1-based for display
        print(f"  e{e+1:<2d}: [{a} {b} {c} {d}]")

    # --- Visualization: temperature field ---
    plt.figure("HW6 Chimney — Temperature Field")
    trisurf_draped_as_contour(nodes, conn, u)
    plt.colorbar(label='Temperature (°C)')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Temperature Distribution (°C) — Q4 FEM')
    plt.tight_layout()

    # --- Node & element numbering overlay ---
    plt.figure("Node & Element Numbering")
    # draw quads as outlines
    for quad in conn:
        poly = nodes[quad, :]
        loop = np.vstack([poly, poly[0]])  # close loop
        plt.plot(loop[:, 0], loop[:, 1], color=(0.4, 0.4, 0.4))
    # node labels
    for n in range(ndof):
        plt.text(nodes[n, 0], nodes[n, 1], f"{n+1}", color='b',
                 ha='center', va='center', fontsize=10)
    # element labels (at centroid)
    for e in range(conn.sha
