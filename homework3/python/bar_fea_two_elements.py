# ========================================================================
# ENGR 4350/6350 — Homework 3
# 1D Bar Finite Element Assembly (Two Linear Elements)
# Trey Brown — Fall 2025
# ========================================================================

import numpy as np

# ----------------------------------------------------------------------
# Problem data
# ----------------------------------------------------------------------
A = 0.5                 # m^2
E = 200e9               # Pa (200 GPa)
b = 10000.0             # N/m (body force)
L_e = 1.0               # element length (m)
AE = A * E

# Two linear elements, nodes at x = 0, 1, 2
n_nodes = 3
n_elems = 2
connectivity = np.array([[1, 2],
                         [2, 3]])  # MATLAB uses 1-based indexing

# ----------------------------------------------------------------------
# Element matrices
# ----------------------------------------------------------------------
k_elem = (AE / L_e) * np.array([[1, -1],
                                [-1, 1]])      # stiffness
f_elem = b * L_e / 2 * np.array([1, 1])        # consistent nodal load

# ----------------------------------------------------------------------
# Global assembly
# ----------------------------------------------------------------------
K = np.zeros((n_nodes, n_nodes))
F = np.zeros(n_nodes)

for e in range(n_elems):
    conn = connectivity[e, :] - 1    # convert to 0-based indexing
    for i_local in range(2):
        i_global = conn[i_local]
        for j_local in range(2):
            j_global = conn[j_local]
            K[i_global, j_global] += k_elem[i_local, j_local]
        F[i_global] += f_elem[i_local]

# ----------------------------------------------------------------------
# Apply Dirichlet BC: u(0) = 0 (node 1 fixed)
# ----------------------------------------------------------------------
fixed_dofs = [0]  # zero-based
free_dofs = [i for i in range(n_nodes) if i not in fixed_dofs]

K_ff = K[np.ix_(free_dofs, free_dofs)]
K_fi = K[np.ix_(free_dofs, fixed_dofs)]
F_f  = F[free_dofs]

u = np.zeros(n_nodes)
u[fixed_dofs] = 0.0

# Reduced system: K_ff * u_f = F_f - K_fi * u_i
RHS = F_f - K_fi @ u[fixed_dofs]
u_free = np.linalg.solve(K_ff, RHS)
u[free_dofs] = u_free

# ----------------------------------------------------------------------
# Compute reactions
# ----------------------------------------------------------------------
reactions = K @ u - F

# ----------------------------------------------------------------------
# Display results
# ----------------------------------------------------------------------
np.set_printoptions(precision=3, suppress=True)
print("\nGlobal stiffness matrix K (N/m):")
print(K)

print("\nGlobal force vector F (N):")
print(F)

print("\nDisplacements u (m) at nodes [x=0, x=1, x=2]:")
print(u)

print("\nReactions R (N) at nodes [x=0, x=1, x=2]:")
print(reactions)

print(f"\nu(x=1) = {u[1]:.6e} m")
print(f"u(x=2) = {u[2]:.6e} m")
