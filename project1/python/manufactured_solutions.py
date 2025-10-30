# manufactured_solution_verification.py
# ------------------------------------------------------------
# Manufactured Solution Verification (1D Tapered Bar)
# ------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Part 1: Manufactured solution and body force
# ------------------------------------------------------------

alpha = 1e-5

# Exact displacement and derivative
def u_exact(x):
    return alpha * x**2

def du_exact(x):
    return 2 * alpha * x

# Bar geometry and material (consistent with tapered_bar_fea)
L = 2.0       # m
E = 210e9     # Pa
A0 = 0.02     # m^2
A1 = 0.005    # m^2

def A_func(x):
    return A0 + (A1 - A0) * (x / L)

# Pack into an input dictionary (similar to MATLAB struct)
inputs = {"L": L, "E": E, "A": A_func}

# Body force function: b(x) = -(1/A) * d/dx(E*A*du/dx)
def body_force(x, inputs):
    dx = 1e-6
    A = inputs["A"]
    term_dA_dx = (A(x + dx) - A(x - dx)) / (2 * dx)
    return -(1.0 / A(x)) * (inputs["E"] * (term_dA_dx * du_exact(x) + A(x) * (2 * alpha)))


# ------------------------------------------------------------
# Part 2: Finite Element Solver with body force
# ------------------------------------------------------------
def fem_solve(inputs, n_el, body_force):
    n_nodes = n_el + 1
    nodes = np.linspace(0, inputs["L"], n_nodes)
    elems = np.array([[i, i + 1] for i in range(n_el)])

    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # Gaussian quadrature (2-point)
    xi_q = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    w_q = [1, 1]

    for e in range(n_el):
        i, j = elems[e]
        x1, x2 = nodes[i], nodes[j]
        le = x2 - x1

        # Element stiffness matrix
        Ae = 0.5 * (inputs["A"](x1) + inputs["A"](x2))
        ke = inputs["E"] * Ae / le * np.array([[1, -1], [-1, 1]])
        K[i:i + 2, i:i + 2] += ke

        # Element load vector (body force)
        fe = np.zeros(2)
        for xi, w in zip(xi_q, w_q):
            # Map xi → x
            x = 0.5 * ((1 - xi) * x1 + (1 + xi) * x2)
            N = np.array([(1 - xi) / 2, (1 + xi) / 2])
            fe += N * inputs["A"](x) * body_force(x, inputs) * (le / 2) * w

        F[i:i + 2] += fe

    # Boundary condition: u(0) = 0
    u = np.zeros(n_nodes)
    fixed_dofs = [0]
    free_dofs = np.arange(1, n_nodes)

    # Reduced system
    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = F[free_dofs]

    # Solve for free displacements
    u_f = np.linalg.solve(K_ff, F_f)
    u[free_dofs] = u_f

    return nodes, u


# ------------------------------------------------------------
# Part 3: Run and compare
# ------------------------------------------------------------

n_el = 10
nodes, u_num = fem_solve(inputs, n_el, body_force)
u_ex = u_exact(nodes)

# Relative L2 error
error = np.linalg.norm(u_num - u_ex) / np.linalg.norm(u_ex)
print(f"Relative L2 error = {error:.6e}")

# ------------------------------------------------------------
# Part 4: Plot comparison
# ------------------------------------------------------------

plt.figure(figsize=(6, 4))
plt.plot(nodes, u_ex * 1e3, "k--", linewidth=1.2, label="Exact (MMS)")
plt.plot(nodes, u_num * 1e3, "ro-", linewidth=1.3, label="FEM", markerfacecolor="r")
plt.xlabel("x (m)")
plt.ylabel("u (mm)")
plt.title("MMS Verification of FEM")
plt.legend(loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.show()
