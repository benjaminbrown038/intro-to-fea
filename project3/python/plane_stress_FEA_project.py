# plane_stress_FEA_project.py
# ENGR 4350/6350 Project 3 — 2D Linear Elasticity (Plane Stress / Plane Strain)
# -------------------------------------------------------------------------
# - 2 DOF per node: ux, uy
# - Q4 (4-node) isoparametric elements
# - Gauss quadrature integration (2x2)
# - Dirichlet BC (fixed displacements)
# - Neumann BC (uniform traction)
# - Post-process: stress tensor, von Mises stress, displacement vectors
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


# ----------------------------
# Material
# ----------------------------
E = 210e9
nu = 0.30
plane_stress = True

if plane_stress:
    D = E / (1 - nu**2) * np.array([
        [1,   nu,  0],
        [nu,  1,   0],
        [0,   0,  (1 - nu) / 2],
    ])
else:
    D = E / ((1 + nu) * (1 - 2 * nu)) * np.array([
        [1 - nu, nu, 0],
        [nu, 1 - nu, 0],
        [0, 0, (1 - 2 * nu) / 2],
    ])


# ----------------------------
# Geometry / Mesh
# ----------------------------
Lx, Ly = 1.0, 0.8
nx, ny = 3, 2
order = 1


@dataclass
class Mesh:
    X: np.ndarray
    Y: np.ndarray
    conn: np.ndarray
    ne: int
    nn: int
    order: int


def make_rect_mesh(Lx, Ly, nx, ny, order=1):
    if order != 1:
        raise ValueError("Only Q4 implemented (set order=1)")
    xs = np.linspace(0, Lx, nx + 1)
    ys = np.linspace(0, Ly, ny + 1)
    Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
    X = Xg.flatten(order="F")
    Y = Yg.flatten(order="F")

    conn = np.zeros((nx * ny, 4), dtype=int)
    e = 0
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i
            n2 = n1 + 1
            n3 = n2 + (nx + 1)
            n4 = n1 + (nx + 1)
            conn[e, :] = [n1, n2, n3, n4]
            e += 1
    return Mesh(X, Y, conn, conn.shape[0], X.size, order)


mesh = make_rect_mesh(Lx, Ly, nx, ny, order)
print(f"Mesh: {mesh.ne} elements, {mesh.nn} nodes")


# ----------------------------
# Boundary Conditions
# ----------------------------
left_nodes = np.where(np.abs(mesh.X) < 1e-12)[0]
fixed_dofs = np.concatenate([[2 * n, 2 * n + 1] for n in left_nodes])
fixed_vals = np.zeros_like(fixed_dofs, dtype=float)

# ----------------------------
# Boundary edges
# ----------------------------
def boundary_edge(mesh, which):
    edge = []
    tol = 1e-12
    Xmax, Ymax = mesh.X.max(), mesh.Y.max()
    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xe, ye = mesh.X[nodes], mesh.Y[nodes]
        if which == "bottom" and np.all(np.abs(ye[[0, 1]] - 0.0) < tol):
            edge.append([e, 1])
        elif which == "right" and np.all(np.abs(xe[[1, 2]] - Xmax) < tol):
            edge.append([e, 2])
        elif which == "top" and np.all(np.abs(ye[[2, 3]] - Ymax) < tol):
            edge.append([e, 3])
        elif which == "left" and np.all(np.abs(xe[[3, 0]] - 0.0) < tol):
            edge.append([e, 4])
    return np.array(edge, dtype=int)


right_edges = boundary_edge(mesh, "right")
traction = lambda x, y: np.array([1e6, 0.0])


# ----------------------------
# Gauss & Shape Functions
# ----------------------------
def gauss1D(n=2):
    if n == 2:
        gp = np.array([-1, 1]) / np.sqrt(3)
        w = np.array([1, 1])
    elif n == 3:
        gp = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
        w = np.array([5/9, 8/9, 5/9])
    else:
        raise ValueError("Use n=2 or 3")
    return gp, w


def shape_Q4():
    gp, w = gauss1D(2)
    xis, etas = np.meshgrid(gp, gp)
    wx, wy = np.meshgrid(w, w)
    wts = (wx * wy).ravel()
    xis = xis.ravel()
    etas = etas.ravel()
    m = len(xis)
    N = np.zeros((m, 4))
    dN_dxi = np.zeros((m, 4))
    dN_deta = np.zeros((m, 4))
    for g in range(m):
        xi, eta = xis[g], etas[g]
        N[g, :] = 0.25 * np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ])
        dN_dxi[g, :] = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta[g, :] = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return N, dN_dxi, dN_deta, wts, xis, etas


# ----------------------------
# Assembly (Elastic)
# ----------------------------
def assemble_elastic(mesh, D, traction, right_edges):
    ndof = 2 * mesh.nn
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    N, dN_dxi, dN_deta, wts, _, _ = shape_Q4()

    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xe, ye = mesh.X[nodes], mesh.Y[nodes]
        Ke = np.zeros((8, 8))
        for g in range(len(wts)):
            J = np.array([
                [dN_dxi[g, :] @ xe, dN_deta[g, :] @ xe],
                [dN_dxi[g, :] @ ye, dN_deta[g, :] @ ye],
            ])
            detJ = np.linalg.det(J)
            if detJ <= 0:
                raise ValueError(f"Non-positive detJ in element {e+1}")
            invJ = np.linalg.inv(J)
            dNdxdy = np.column_stack([dN_dxi[g, :], dN_deta[g, :]]) @ invJ
            B = np.zeros((3, 8))
            for i in range(4):
                B[:, 2*i:2*i+2] = np.array([
                    [dNdxdy[i, 0], 0],
                    [0, dNdxdy[i, 1]],
                    [dNdxdy[i, 1], dNdxdy[i, 0]]
                ])
            Ke += (B.T @ D @ B) * detJ * wts[g]
        dofs = np.ravel(np.column_stack([2*nodes, 2*nodes + 1]))
        K[np.ix_(dofs, dofs)] += Ke

    # Traction on right edges
    for e, ledge in right_edges:
        nodes2 = mesh.conn[e, [1, 2]]
        xe, ye = mesh.X[nodes2], mesh.Y[nodes2]
        gp, w = gauss1D(2)
        Ledge = np.hypot(xe[1] - xe[0], ye[1] - ye[0])
        for s, wt in zip(gp, w):
            N1, N2 = (1 - s) / 2, (1 + s) / 2
            xg = N1 * xe[0] + N2 * xe[1]
            yg = N1 * ye[0] + N2 * ye[1]
            t = traction(xg, yg)
            Fe_local = np.vstack([N1 * t, N2 * t]) * (Ledge / 2 * wt)
            dofs = [2*nodes2[0], 2*nodes2[0]+1, 2*nodes2[1], 2*nodes2[1]+1]
            F[dofs] += Fe_local.ravel()
    return K, F


# ----------------------------
# Assembly and Solve
# ----------------------------
K, F = assemble_elastic(mesh, D, traction, right_edges)

ndof = 2 * mesh.nn
free = np.setdiff1d(np.arange(ndof), fixed_dofs)
u = np.zeros(ndof)
F_mod = F[free] - K[np.ix_(free, fixed_dofs)] @ fixed_vals
K_mod = K[np.ix_(free, free)]
u[free] = np.linalg.solve(K_mod, F_mod)


# ----------------------------
# Post Processing
# ----------------------------
def stress_post(mesh, D, u):
    N, dN_dxi, dN_deta, wts, xis, etas = shape_Q4()
    stress, vonMises, coord_gp = [], [], []
    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xe, ye = mesh.X[nodes], mesh.Y[nodes]
        ue = u[np.ravel(np.column_stack([2*nodes, 2*nodes+1]))]
        for g in range(len(wts)):
            J = np.array([
                [dN_dxi[g, :] @ xe, dN_deta[g, :] @ xe],
                [dN_dxi[g, :] @ ye, dN_deta[g, :] @ ye],
            ])
            invJ = np.linalg.inv(J)
            dNdxdy = np.column_stack([dN_dxi[g, :], dN_deta[g, :]]) @ invJ
            B = np.zeros((3, 8))
            for i in range(4):
                B[:, 2*i:2*i+2] = np.array([
                    [dNdxdy[i, 0], 0],
                    [0, dNdxdy[i, 1]],
                    [dNdxdy[i, 1], dNdxdy[i, 0]]
                ])
            sig = D @ (B @ ue)
            vm = np.sqrt(sig[0]**2 - sig[0]*sig[1] + sig[1]**2 + 3*sig[2]**2)
            xi, eta = xis[g], etas[g]
            Ng = 0.25 * np.array([
                (1 - xi)*(1 - eta),
                (1 + xi)*(1 - eta),
                (1 + xi)*(1 + eta),
                (1 - xi)*(1 + eta)
            ])
            xg = Ng @ xe
            yg = Ng @ ye
            coord_gp.append([xg, yg])
            stress.append(sig)
            vonMises.append(vm)
    return np.array(stress), np.array(vonMises), np.array(coord_gp)


stress, vonMises, coord_gp = stress_post(mesh, D, u)


# ----------------------------
# Plot Results
# ----------------------------
def plot_results(mesh, u, coord_gp, stress, vonMises):
    ux, uy = u[0::2], u[1::2]

    # Node numbering
    plt.figure("Node numbering")
    plt.scatter(mesh.X, mesh.Y, s=25)
    for i in range(mesh.nn):
        plt.text(mesh.X[i], mesh.Y[i], str(i + 1), fontsize=8)
    plt.axis("equal")
    plt.title("Global node numbers")
    plt.xlabel("x"); plt.ylabel("y"); plt.grid(True)

    # Displacement vectors
    umax = np.max(np.sqrt(ux**2 + uy**2)) + 1e-12
    scale = 0.1 * (mesh.X.max() - mesh.X.min()) / umax
    plt.figure("Displacement vectors")
    plt.quiver(mesh.X, mesh.Y, ux, uy, scale=scale)
    plt.axis("equal"); plt.title("Displacement vectors (scaled)")
    plt.xlabel("x"); plt.ylabel("y"); plt.grid(True)

    # Von Mises contour
    from scipy.interpolate import griddata
    xg = np.linspace(mesh.X.min(), mesh.X.max(), 100)
    yg = np.linspace(mesh.Y.min(), mesh.Y.max(), 100)
    Xg, Yg = np.meshgrid(xg, yg)
    Z = griddata(coord_gp, vonMises, (Xg, Yg), method="linear")
    plt.figure("Von Mises stress")
    plt.contourf(Xg, Yg, Z, 20, cmap="jet")
    plt.colorbar(label="Von Mises (Pa)")
    plt.axis("equal"); plt.title("Von Mises stress field")
    plt.xlabel("x [m]"); plt.ylabel("y [m]"); plt.grid(True)
    plt.show()


plot_results(mesh, u, coord_gp, stress, vonMises)


# ----------------------------
# Displacement table & CSV
# ----------------------------
Ux = u[0::2]
Uy = u[1::2]
Utab = np.column_stack([mesh.X, mesh.Y, Ux, Uy])
print("\nDisplacements at each node [ux, uy] (m):")
print("   x          y          ux           uy")
for i in range(mesh.nn):
    print(f"{mesh.X[i]:8.4f} {mesh.Y[i]:8.4f} {Ux[i]:12.6e} {Uy[i]:12.6e}")

np.savetxt("nodal_displacements.csv", Utab,
           delimiter=",", header="x,y,ux,uy", comments="")
print("Saved nodal displacements to nodal_displacements.csv")
