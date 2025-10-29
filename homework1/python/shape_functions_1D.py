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
plt.subplot(1, 2, 1)
plt.plot(xi, u_lin, 'b', linewidth=1.5)
plt.plot([-1, 1], d_lin, 'ro', markerfacecolor='r')
plt.title('u(ξ) = ξ + 1')
plt.xlabel('ξ')
plt.ylabel('u(ξ)')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(xi, du_lin, 'r', linewidth=1.5)
plt.title('du/dξ = 1')
plt.xlabel('ξ')
plt.ylabel('du/dξ')
plt.grid(True)

plt.tight_layout()
plt.show()

# ======================================================================
# 3-node Quadratic Element
# ======================================================================
d_quad = np.array([-1, -1, 1])   # nodal values at xi = -1, 0, 1

u_quad = np.zeros_like(xi)
du_quad = np.zeros_like(xi)

for k in range(len(xi)):
    N = np.array([
        0.5 * xi[k] * (xi[k] - 1),
        1 - xi[k]**2,
        0.5 * xi[k] * (xi[k] + 1)
    ])
    B = np.array([
        xi[k] - 0.5,
        -2 * xi[k],
        xi[k] + 0.5
    ])
    u_quad[k] = N @ d_quad
    du_quad[k] = B @ d_quad

plt.figure(figsiz
