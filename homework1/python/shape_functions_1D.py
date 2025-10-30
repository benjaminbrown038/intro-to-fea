# ========================================================================
# ENGR 4350/6350 — Homework 1
# Finite Element Shape Functions Validation (1D)
# Trey Brown — Fall 2025
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1) Setup
# ----------------------------------------------------------------------

xi = np.linspace(-1, 1, 200)

# ======================================================================
# 2-node Linear Element
# ======================================================================
d_lin = np.array([0, 2])   # nodal values at xi = -1 and xi = 1

u_lin = np.zeros_like(xi)
du_lin = np.zeros_like(xi)

for k in range(len(xi)):
    N = np.array([(1 - xi[k]) / 2, (1 + xi[k]) / 2])
    B = np.array([-0.5, 0.5])
    u_lin[k] = N @ d_lin
    du_lin[k] = B @ d_lin

plt.figure(figsize=(10, 4))

# Left subplot: Linear shape function
plt.subplot(1, 2, 1)
plt.plot(xi, u_lin, 'b', linewidth=1.5)
plt.plot([-1, 1], d_lin, 'ro', markerfacecolor='r')
plt.title('Linear Element Shape Function: $u(\\xi)$')
plt.xlabel('Parent Coordinate $\\xi$')
plt.ylabel('Displacement $u(\\xi)$')
plt.grid(True)

# Right subplot: Derivative
plt.subplot(1, 2, 2)
plt.plot(xi, du_lin, 'r', linewidth=1.5)
plt.title('Derivative of Shape Function: $du/d\\xi$')
plt.xlabel('Parent Coordinate $\\xi$')
plt.ylabel('Slope $du/d\\xi$')
plt.grid(True)

plt.suptitle('Two-Node Linear Finite Element — Shape Function and Derivative', fontsize=13)
plt.tight_layout()
plt.show()


# ======================================================================
# 3-node Quadratic Element
# ======================================================================
d_quad = np.array([-1, -1, 1])   # nodal values at xi = -1, 0, 1

u_quad = np.zeros_like(xi)
du_quad = np.zeros_like(xi)

for k in range(len(xi)):
    # Shape functions
    N = np.array([
        0.5 * xi[k] * (xi[k] - 1),
        1 - xi[k]**2,
        0.5 * xi[k] * (xi[k] + 1)
    ])
    # Derivatives
    B = np.array([
        xi[k] - 0.5,
        -2 * xi[k],
        xi[k] + 0.5
    ])
    u_quad[k] = N @ d_quad
    du_quad[k] = B @ d_quad

# ======================================================================
# Plot results for the 3-node element
# ======================================================================
plt.figure(figsize=(10, 4))

# Left subplot: Quadratic shape function
plt.subplot(1, 2, 1)
plt.plot(xi, u_quad, 'b', linewidth=1.5)
plt.plot([-1, 0, 1], d_quad, 'ro', markerfacecolor='r')
plt.title('Quadratic Element Shape Function')
plt.xlabel('Parent Coordinate $\\xi$ (Natural Domain)')
plt.ylabel('Interpolated Displacement $u(\\xi)$')
plt.grid(True)

# Right subplot: Derivative
plt.subplot(1, 2, 2)
plt.plot(xi, du_quad, 'r', linewidth=1.5)
plt.title('Derivative of Quadratic Shape Function')
plt.xlabel('Parent Coordinate $\\xi$ (Natural Domain)')
plt.ylabel('Shape Function Gradient $du/d\\xi$')
plt.grid(True)

plt.suptitle('Three-Node Quadratic Finite Element — Shape Function and Derivative', fontsize=13)
plt.tight_layout()
plt.show()


# ======================================================================
#  End of file
# ======================================================================
