# tapered_bar_fea.py
# ------------------------------------------------------------
# 1D Finite Element Analysis of a Tapered Bar
# ------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Part 1: Inputs
# ----------------------------

# Bar properties
L = 2.0                # Length (m)
E = 210e9              # Young's modulus (Pa)
A0 = 0.02              # Area at x = 0 (m^2)
A1 = 0.005             # Area at x = L (m^2)

u0_prescribed = 0.0    # Fixed end displacement (m)
sigma_L_prescribed = None  # No stress specified at x = L

# Area variation functions
def linear_taper(A0, A1, L):
    return lambda x: A0 + (A1 - A0) * (x / L)

def exponential_taper(A0, A1, L):
    return lambda x: A0 * (A1 / A0) ** (x / L)

# Choose area function
A_func = linear_taper(A0, A1, L)

# --- Summary output ---
print("Tapered bar inputs summary")
print("--------------------------")
print(f"L = {L:.3f} m")
print(f"E = {E:.3e} Pa")
print(f"u(0) = {u0_prescribed:.3f}")
if sigma_L_prescribed is None:
    print("sigma(L) = None")
else:
    print(f"sigma(L) = {sigma_L_prescribed:.3e}")

xs = np.linspace(0, L, 5)
areas = A_func(xs)
print("Sample A(x):")
for x, A in zip(xs, areas):
    print(f"  A({x:.4g}) = {A:.6g}")

# ----------------------------
# Part 2: Finite Element Solver
# ----------------------------

n_el = 10
n_nodes = n_el + 1
nodes = np.linspace(0, L, n_nodes)
elems = np.array([[i, i + 1] for i in range(n_el)])

K = np.zeros((n_nodes, n_nodes))
F = np.zeros(n_nodes)

# Assemble global stiffness matrix
for e in range(n_el):
    i, j = elems[e]
    x1, x2 = nodes[i], nodes[j]
    le = x2 - x1
    Ae = (A_func(x1) + A_func(x2)) / 2
    ke = E * Ae / le * np.array([[1, -1], [-1, 1]])
    K[i:i + 2, i:i + 2] += ke

# Apply boundary conditions
fixed_dofs = [0]
free_dofs = np.arange(1, n_nodes)

F[-1] = 1e5   # Force at the right end (N)

# Solve
K_ff = K[np.ix_(free_dofs, free_dofs)]
F_f = F[free_dofs]
u = np.zeros(n_nodes)
u_f = np.linalg.solve(K_ff, F_f)
u[free_dofs] = u_f

# ----------------------------
# Part 3: Stresses
# ----------------------------

stresses = np.zeros(n_el)
for e in range(n_el):
    i, j = elems[e]
    le = nodes[j] - nodes[i]
    du = u[j] - u[i]
    strain = du / le
    stresses[e] = E * strain

# ----------------------------
# Part 4: Plots
# ----------------------------

fig, axs = plt.subplots(1, 3, figsize=(14, 4))

# (a) Area distribution
xs_dense = np.linspace(0, L, 200)
axs[0].plot(xs_dense, A_func(xs_dense), linewidth=1.5)
axs[0].set_title("Cross-sectional Area A(x)")
axs[0].set_xlabel("x (m)")
axs[0].set_ylabel("A(x) (m²)")
axs[0].grid(True)

# (b) Displacements
axs[1].plot(nodes, u * 1e3, 'o-', linewidth=1.5)
axs[1].set_title("Nodal Displacements")
axs[1].set_xlabel("x (m)")
axs[1].set_ylabel("u (mm)")
axs[1].grid(True)

# (c) Stresses
x_centers = nodes[:-1] + np.diff(nodes) / 2
axs[2].plot(x_centers, stresses * 1e-6, 's-r', linewidth=1.5)
axs[2].set_title("Element Stresses")
axs[2].set_xlabel("x (m)")
axs[2].set_ylabel("Stress (MPa)")
axs[2].grid(True)

plt.suptitle("Tapered Bar Finite Element Results")
plt.tight_layout()
plt.show()
