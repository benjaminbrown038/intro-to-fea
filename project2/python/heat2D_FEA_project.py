# heat2D_FEA_project.py
# ENGR 4350/6350 — Project 2: 2-D steady heat conduction (FEM)
# 1) Inputs
# 2) Dirichlet BCs
# 3) Neumann BCs
# 4) Mesh
# 5) Assembly
# 6) Boundary Conditions
# 7) Post Processing
# 8) Mesh Generation
# 9) Element / Shape Function
# 10) Assembly (element → global)
# 11) Boundary Helpers
# 12) Dirichlet Solver

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# ----------------------------
# 1) Inputs
# ----------------------------
k   = 20.0
Lx  = 1.0
Ly  = 0.8
order = 1       # 1 = Q4, 2 = Q9
nx, ny = 3, 2   # elements in x, y

# Dirichlet BC (left)
T_left    = 400.0
use_Tleft = True

# Neumann BC (top edge)
q_top     = 1.0e4
use_qtop  = True


# ----------------------------
# Mesh data structure
# ----------------------------
@dataclass
class Mesh:
    X: np.ndarray           # (nn,)
    Y: np.ndarray           # (nn,)
    conn: np.ndarray        # (ne, nper) with nper=4 (Q4) or 9 (Q9)
    nn: int
    ne: int
    order: int
    nx: int
    ny: int
    nxn: int                # number of nodes in x
    nyn: int                # number of nodes in y


# ----------------------------
# 8) Mesh Generation
# ----------------------------
def make_rect_mesh(Lx: float, Ly: float, nx: int, ny: int, order: int) -> Mesh:
    if order == 1:
        xs = np.linspace(0.0, Lx, nx + 1)
        ys = np.linspace(0.0, Ly, ny + 1)
        nxn, nyn = nx + 1, ny + 1
        # Mimic MATLAB ndgrid indexing
        Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
        # Flatten like MATLAB (:), i.e., column-major
        X = Xg.flatten(order="F")
        Y = Yg.flatten(order="F")

        def node_id(i, j):  # zero-based i in [0,nxn-1], j in [0,nyn-1]
            return j * nxn + i

        conn = np.zeros((nx * ny, 4), dtype=int)
        e = 0
        for j in range(ny):
            for i in range(nx):
                n1 = node_id(i,     j)
                n2 = node_id(i + 1, j)
                n3 = node_id(i + 1, j + 1)
                n4 = node_id(i,     j + 1)
                conn[e, :] = [n1, n2, n3, n4]
                e += 1

    elif order == 2:
        xs = np.linspace(0.0, Lx, 2 * nx + 1)
        ys = np.linspace(0.0, Ly, 2 * ny + 1)
        nxn, nyn = 2 * nx + 1, 2 * ny + 1
        Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
        X = Xg.flatten(order="F")
        Y = Yg.flatten(order="F")

        def node_id(i, j):
            return j * nxn + i

        conn = np.zeros((nx * ny, 9), dtype=int)
        e = 0
        for j in range(ny):
            for i in range(nx):
                i0 = 2 * i
                j0 = 2 * j
                # Corner + midside + center with the same pattern as MATLAB
                n11 = node_id(i0,   j0)
                n12 = node_id(i0+1, j0)
                n13 = node_id(i0+2, j0)
                n21 = node_id(i0,   j0+1)
                n22 = node_id(i0+1, j0+1)
                n23 = node_id(i0+2, j0+1)
                n31 = node_id(i0,   j0+2)
                n32 = node_id(i0+1, j0+2)
                n33 = node_id(i0+2, j0+2)
                # MATLAB ordering: [n11 n13 n33 n31 n12 n23 n32 n21 n22]
                conn[e, :] = [n11, n13, n33, n31, n12, n23, n32, n21, n22]
                e += 1
    else:
        raise ValueError("order must be 1 (Q4) or 2 (Q9)")

    nn = X.size
    ne = conn.shape[0]
    return Mesh(X=X, Y=Y, conn=conn, nn=nn, ne=ne, order=order, nx=nx, ny=ny, nxn=nxn, nyn=nyn)


# ----------------------------
# 9) Element / Shape Function
# ----------------------------
def gauss1D(n: int):
    if n == 1:
        return np.array([0.0]), np.array([2.0])
    elif n == 2:
        v = 1.0 / np.sqrt(3.0)
        return np.array([-v, v]), np.array([1.0, 1.0])
    elif n == 3:
        r = np.sqrt(3.0 / 5.0)
        return np.array([-r, 0.0, r]), np.array([5/9, 8/9, 5/9])
    else:
        raise ValueError("use n=1..3")


def element_shape(order: int, rule: int | None = None):
    if order == 1:
        if rule is None:
            rule = 2
        gp1, w1 = gauss1D(rule)
        xis, etas = np.meshgrid(gp1, gp1, indexing="xy")
        wx, wy = np.meshgrid(w1, w1, indexing="xy")
        wts = (wx * wy).ravel()
        xis = xis.ravel()
        etas = etas.ravel()
        m = xis.size
        N = np.zeros((m, 4))
        dN_dxi = np.zeros((m, 4))
        dN_deta = np.zeros((m, 4))
        for g in range(m):
            xi = xis[g]
            eta = etas[g]
            N[g, :] = 0.25 * np.array([
                (1 - xi) * (1 - eta),
                (1 + xi) * (1 - eta),
                (1 + xi) * (1 + eta),
                (1 - xi) * (1 + eta),
            ])
            dN_dxi[g, :]  = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
            dN_deta[g, :] = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi),  (1-xi)])
        return N, dN_dxi, dN_deta, wts, xis, etas

    elif order == 2:
        if rule is None:
            rule = 3
        gp1, w1 = gauss1D(rule)
        xis, etas = np.meshgrid(gp1, gp1, indexing="xy")
        wx, wy = np.meshgrid(w1, w1, indexing="xy")
        wts = (wx * wy).ravel()
        xis = xis.ravel()
        etas = etas.ravel()
        m = xis.size
        N = np.zeros((m, 9))
        dN_dxi = np.zeros((m, 9))
        dN_deta = np.zeros((m, 9))

        def L(t):
            return np.array([0.5 * t * (t - 1.0), 1.0 - t**2, 0.5 * t * (t + 1.0)])

        def dL(t):
            return np.array([t - 0.5, -2.0 * t, t + 0.5])

        # map [r,c] -> id like MATLAB M = [[1 5 2],[8 9 6],[4 7 3]]
        M = np.array([[1, 5, 2],
                      [8, 9, 6],
                      [4, 7, 3]], dtype=int) - 1  # zero-based

        for g in range(m):
            xi = xis[g]
            eta = etas[g]
            Lx, Ly = L(xi), L(eta)
            dLx, dLy = dL(xi), dL(eta)
            for r in range(3):
                for c in range(3):
                    idx = M[r, c]
                    N[g, idx]       = Lx[c] * Ly[r]
                    dN_dxi[g, idx]  = dLx[c] * Ly[r]
                    dN_deta[g, idx] = Lx[c] * dLy[r]
        return N, dN_dxi, dN_deta, wts, xis, etas

    else:
        raise ValueError("order must be 1 (Q4) or 2 (Q9)")


# ----------------------------
# 10) Assembly (element → global)
# ----------------------------
def assemble_global_conduction(mesh: Mesh, k: float):
    nn, ne = mesh.nn, mesh.ne
    K = np.zeros((nn, nn))
    F = np.zeros(nn)

    rule = 2 if mesh.order == 1 else 3
    N, dN_dxi, dN_deta, wts, _, _ = element_shape(mesh.order, rule)
    m = wts.size

    for e in range(ne):
        nodes = mesh.conn[e, :]
        xe = mesh.X[nodes]
        ye = mesh.Y[nodes]
        Ke = np.zeros((nodes.size, nodes.size))
        for g in range(m):
            J = np.array([
                [dN_dxi[g, :].dot(xe),  dN_deta[g, :].dot(xe)],
                [dN_dxi[g, :].dot(ye),  dN_deta[g, :].dot(ye)],
            ])
            detJ = np.linalg.det(J)
            if detJ <= 0:
                raise ValueError(f"Non-positive detJ in element {e+1}")
            invJ = np.linalg.inv(J)
            dNdxdy = np.column_stack([dN_dxi[g, :], dN_deta[g, :]]).dot(invJ)  # (nper,2)
            B = dNdxdy.T  # (2, nper)
            Ke += (k * (B.T @ B)) * detJ * wts[g]

        # Assemble
        K[np.ix_(nodes, nodes)] += Ke

    return K, F


# ----------------------------
# 11) Boundary Helpers
# ----------------------------
def boundary_edge(mesh: Mesh, which: str) -> np.ndarray:
    edge = []
    tol = 1e-12
    Xmax = mesh.X.max()
    Ymax = mesh.Y.max()
    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xe = mesh.X[nodes]
        ye = mesh.Y[nodes]
        if which == "bottom":
            if np.all(np.abs(ye[[0, 1]] - 0.0) < tol):
                edge.append([e, 1])
        elif which == "right":
            if np.all(np.abs(xe[[1, 2]] - Xmax) < tol):
                edge.append([e, 2])
        elif which == "top":
            if np.all(np.abs(ye[[2, 3]] - Ymax) < tol):
                edge.append([e, 3])
        elif which == "left":
            if np.all(np.abs(xe[[3, 0]] - 0.0) < tol):
                edge.append([e, 4])
        else:
            raise ValueError("which ∈ {bottom,right,top,left}")
    return np.array(edge, dtype=int)


def apply_edge_flux(F: np.ndarray, mesh: Mesh, edge_list: np.ndarray, qn_fun):
    gp, wts = gauss1D(2)
    for r in range(edge_list.shape[0]):
        e, ledge = int(edge_list[r, 0]), int(edge_list[r, 1])
        nodes = mesh.conn[e, :]
        xe = mesh.X[nodes]
        ye = mesh.Y[nodes]

        for g in range(gp.size):
            s = gp[g]
            wt = wts[g]
            if ledge == 1:   # bottom: (1-2)
                N = np.array([(1 - s) / 2, (1 + s) / 2, 0, 0])
                ids = [0, 1]
            elif ledge == 2: # right: (2-3)
                N = np.array([0, (1 - s) / 2, (1 + s) / 2, 0])
                ids = [1, 2]
            elif ledge == 3: # top: (3-4)
                N = np.array([0, 0, (1 - s) / 2, (1 + s) / 2])
                ids = [2, 3]
            elif ledge == 4: # left: (4-1)
                N = np.array([(1 + s) / 2, 0, 0, (1 - s) / 2])
                ids = [3, 0]
            else:
                raise ValueError("ledge must be 1..4")

            xg = float(N @ xe[:4])
            yg = float(N @ ye[:4])
            qn = float(qn_fun(xg, yg))
            dxds = 0.5 * (xe[ids[1]] - xe[ids[0]])
            dyds = 0.5 * (ye[ids[1]] - ye[ids[0]])
            Jedge = np.hypot(dxds, dyds)
            F[nodes[:4]] += N * (qn * Jedge * wt)
    return F


# ----------------------------
# 12) Dirichlet Solver
# ----------------------------
def solve_with_dirichlet(K: np.ndarray, F: np.ndarray, bc: dict):
    nn = K.shape[0]
    T = np.full(nn, np.nan)
    fixed = np.unique(np.asarray(bc.get("dir_nodes", []), dtype=int))
    if fixed.size:
        T[fixed] = np.asarray(bc.get("dir_vals", []), dtype=float)
    free = np.setdiff1d(np.arange(nn), fixed)

    F_mod = F[free] - K[np.ix_(free, fixed)] @ T[fixed] if fixed.size else F[free]
    K_mod = K[np.ix_(free, free)]
    T[free] = np.linalg.solve(K_mod, F_mod)

    fixed_pack = {"fixed": fixed, "Tfixed": T[fixed]}
    return T, free, fixed_pack


# ----------------------------
# Post-processing helpers
# ----------------------------
def post_flux_on_gauss(mesh: Mesh, k: float, T: np.ndarray):
    rule = 2 if mesh.order == 1 else 3
    _, dN_dxi, dN_deta, wts, xis, etas = element_shape(mesh.order, rule)
    gpXY = []
    fluxXY = []
    dTdx_gp = []
    dTdy_gp = []

    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xe = mesh.X[nodes]
        ye = mesh.Y[nodes]
        Te = T[nodes]
        for g in range(wts.size):
            J = np.array([
                [dN_dxi[g, :].dot(xe),  dN_deta[g, :].dot(xe)],
                [dN_dxi[g, :].dot(ye),  dN_deta[g, :].dot(ye)],
            ])
            invJ = np.linalg.inv(J)
            dNdxdy = np.column_stack([dN_dxi[g, :], dN_deta[g, :]]).dot(invJ)  # (nper,2)
            gradT = dNdxdy.T @ Te  # (2,)
            q = -k * gradT

            xi, eta = xis[g], etas[g]
            # Q4 interpolation of the Gauss-point location (corner N)
            N4 = 0.25 * np.array([
                (1 - xi) * (1 - eta),
                (1 + xi) * (1 - eta),
                (1 + xi) * (1 + eta),
                (1 - xi) * (1 + eta),
            ])
            # Use first 4 (corner) nodes for location
            xg = float(N4 @ xe[:4])
            yg = float(N4 @ ye[:4])

            gpXY.append([xg, yg])
            fluxXY.append([q[0], q[1]])
            dTdx_gp.append(gradT[0])
            dTdy_gp.append(gradT[1])

    return (np.array(gpXY), np.array(fluxXY),
            np.array(dTdx_gp), np.array(dTdy_gp))


def report_and_plots(mesh: Mesh, T: np.ndarray,
                     gpXY: np.ndarray, fluxXY: np.ndarray,
                     dTdx_gp: np.ndarray, dTdy_gp: np.ndarray,
                     order: int):

    # Numbering plot
    plt.figure("Numbering")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.scatter(mesh.X, mesh.Y, s=25)
    for n in range(mesh.nn):
        plt.text(mesh.X[n], mesh.Y[n], str(n + 1), fontsize=8)
    for e in range(mesh.ne):
        nodes = mesh.conn[e, :]
        xc = float(np.mean(mesh.X[nodes]))
        yc = float(np.mean(mesh.Y[nodes]))
        plt.text(xc, yc, f"e{e+1}", color=(0.85, 0.0, 0.2), fontweight="bold")
    plt.title(f"Global node numbers and element labels ({'Q4' if order==1 else 'Q9'})")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.grid(True)

    # Temperature contour (structured grid → pcolormesh)
    plt.figure("Temperature Contour")
    # Reshape with MATLAB-like memory order
    Xg = mesh.X.reshape((mesh.nxn, mesh.nyn), order="F")
    Yg = mesh.Y.reshape((mesh.nxn, mesh.nyn), order="F")
    Tg = T.reshape((mesh.nxn, mesh.nyn), order="F")
    # pcolormesh expects cell edges; use node grid directly for smooth look
    pcm = plt.pcolormesh(Xg, Yg, Tg, shading="gouraud")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.colorbar(pcm, label="T")
    plt.title("Temperature field T(x,y)")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.grid(True, alpha=0.3)

    # Heat flux quiver at Gauss points
    plt.figure("Heat Flux")
    plt.quiver(gpXY[:, 0], gpXY[:, 1], fluxXY[:, 0], fluxXY[:, 1], angles="xy")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title("Heat-flux vectors q = -k ∇T (Gauss points)")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.grid(True)

    plt.tight_layout()
    plt.show()


def disp_table_T(mesh: Mesh, T: np.ndarray):
    print("\n=== Nodal Temperatures (K) ===")
    header = f"{'Node':>6} {'x':>12} {'y':>12} {'T':>14}"
    print(header)
    for i in range(mesh.nn):
        print(f"{i+1:6d} {mesh.X[i]:12.6g} {mesh.Y[i]:12.6g} {T[i]:14.6g}")


# ----------------------------
# Main driver (equivalent to MATLAB top-level function)
# ----------------------------
def main():
    mesh = make_rect_mesh(Lx, Ly, nx, ny, order)
    print(f"Mesh: {mesh.ne} elements, {mesh.nn} nodes (order={'Q4' if order==1 else 'Q9'})")

    # Assembly
    K, F = assemble_global_conduction(mesh, k)

    # Boundary Conditions (Dirichlet left)
    bc = {}
    if use_Tleft:
        left_nodes = np.where(np.abs(mesh.X) < 1e-12)[0]
        bc["dir_nodes"] = left_nodes
        bc["dir_vals"]  = T_left * np.ones(left_nodes.size)
    else:
        bc["dir_nodes"] = np.array([], dtype=int)
        bc["dir_vals"]  = np.array([], dtype=float)

    # Neumann BC on top
    if use_qtop:
        top_edge = boundary_edge(mesh, "top")
        F = apply_edge_flux(F, mesh, top_edge, lambda x, y: q_top)

    # Solve with Dirichlet
    T, free, fixed_pack = solve_with_dirichlet(K, F, bc)

    # Post-processing
    gpXY, fluxXY, dTdx_gp, dTdy_gp = post_flux_on_gauss(mesh, k, T)
    report_and_plots(mesh, T, gpXY, fluxXY, dTdx_gp, dTdy_gp, order)
    disp_table_T(mesh, T)

    # Save nodal temperatures
    out = np.column_stack([np.arange(1, mesh.nn + 1), mesh.X, mesh.Y, T])
    np.savetxt("nodal_temperatures.csv", out,
               fmt=["%d", "%.10g", "%.10g", "%.10g"],
               delimiter=",",
               header="Node,x,y,T", comments="")
    print("Saved nodal temperatures to nodal_temperatures.csv")


if __name__ == "__main__":
    main()
