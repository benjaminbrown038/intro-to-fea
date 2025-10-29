# ========================================================================
# ENGR 4350/6350 — Homework 5
# Compute element stiffness [Ke] and heat generation flux {fe}
# for one 4-node quadrilateral element using 2x2 Gauss quadrature
# Trey Brown — Fall 2025
# ========================================================================

import numpy as np

# ----------------------------------------------------------------------
# Given parameters
# ----------------------------------------------------------------------
k = 5.0        # thermal conductivity
s = 10.0       # heat generation
xy = np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [0.5, 0.8],
    [0.0, 0.5]
])  # nodal coordinates [x_i, y_i]

# ----------------------------------------------------------------------
# 2×2 Gauss quadrature setup
# ----------------------------------------------------------------------
gp = np.array([-1, 1]) / np.sqrt(3)
w  = np.array([1.0, 1.0])

# Initialize results
Ke = np.zeros((4, 4))
fe = np.zeros(4)

# ----------------------------------------------------------------------
# Loop over Gauss points
# ----------------------------------------------------------------------
for i in range(2):
    for j in range(2):
        xi  = gp[i]
        eta = gp[j]
        wi  = w[i]
        wj  = w[j]

        # --- Shape functions (Q4) ---
        N = 0.25 * np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ])

        # --- Derivatives wrt parent coords ---
        dN_dxi  = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi),  (1 - xi)])

        # --- Jacobian matrix ---
        J = np.vstack((dN_dxi, dN_deta)) @ xy
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        # --- Derivatives wrt global coords ---
        dN = invJ @ np.vstack((dN_dxi, dN_deta))
        dNdx = dN[0, :]
        dNdy = dN[1, :]

        # --- B matrix (2×4) ---
        B = np.vstack((dNdx, dNdy))

        # --- Element stiffness contribution ---
        Ke += (B.T @ (k * B)) * detJ * wi * wj

        # --- Flux vector due to heat generation ---
        fe += (N * s) * detJ * wi * wj

# ----------------------------------------------------------------------
# Display results
# ----------------------------------------------------------------------
np.set_printoptions(precision=6, suppress=True)

print("Element stiffness matrix [Ke] (4×4):")
print(Ke)

print("\nFlux vector due to heat generation {fe}:")
print(fe)
