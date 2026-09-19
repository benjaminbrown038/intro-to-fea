# Tier 1 – Basic FEA Foundations

## Overview

Analyze a 2D cantilever beam under a uniformly distributed
traction on its free end.

The simulation computes nodal displacements, element strains,
element stresses, and the deformed shape.

## Physical Assumptions

- Static equilibrium.
- Small strains and small rotations.
- Linear elastic, homogeneous, isotropic material.
- Plane stress with constant out-of-plane thickness.
- A fully fixed left edge.
- A uniformly loaded right edge.

## Geometry and Displacement

The rectangular domain is:

$$
\Omega=[0,L]\times[0,H].
$$

Here, $L$ is beam length, $H$ is beam height, and $t$ is
out-of-plane thickness.

The displacement field is:

$$
\mathbf{u}(x,y)
=
\begin{bmatrix}
u_x(x,y)\\
u_y(x,y)
\end{bmatrix}.
$$

## Governing Equations

### Static Equilibrium

With no body forces:

$$
\nabla\cdot\boldsymbol{\sigma}=\mathbf{0}.
$$

In component form:

$$
\frac{\partial\sigma_{xx}}{\partial x}
+
\frac{\partial\tau_{xy}}{\partial y}
=0,
$$

$$
\frac{\partial\tau_{xy}}{\partial x}
+
\frac{\partial\sigma_{yy}}{\partial y}
=0.
$$

### Strain–Displacement Relationship

Using engineering shear strain:

$$
\boldsymbol{\varepsilon}_v
=
\begin{bmatrix}
\varepsilon_{xx}\\
\varepsilon_{yy}\\
\gamma_{xy}
\end{bmatrix}
=
\begin{bmatrix}
\frac{\partial u_x}{\partial x}\\
\frac{\partial u_y}{\partial y}\\
\frac{\partial u_x}{\partial y}
+
\frac{\partial u_y}{\partial x}
\end{bmatrix}.
$$

Engineering shear strain satisfies
$\gamma_{xy}=2\varepsilon_{xy}$.

### Plane-Stress Constitutive Relationship

$$
\boldsymbol{\sigma}_v
=
\mathbf{D}\boldsymbol{\varepsilon}_v,
\qquad
\boldsymbol{\sigma}_v
=
\begin{bmatrix}
\sigma_{xx}\\
\sigma_{yy}\\
\tau_{xy}
\end{bmatrix}.
$$

$$
\mathbf{D}
=
\frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0\\
\nu & 1 & 0\\
0 & 0 & \frac{1-\nu}{2}
\end{bmatrix}.
$$

Here, $E$ is Young's modulus and $\nu$ is Poisson's ratio.
Plane stress assumes:

$$
\sigma_{zz}=\tau_{xz}=\tau_{yz}=0.
$$

## Boundary Conditions

### Fixed Edge

At $x=0$:

$$
u_x=0,
\qquad
u_y=0.
$$

### Uniform End Traction

At $x=L$, apply a downward traction of magnitude $q$:

$$
\boldsymbol{\sigma}\mathbf{n}
=
\begin{bmatrix}
0\\
-q
\end{bmatrix}.
$$

Here, $q$ is force per unit area. The total applied force
magnitude is:

$$
P=qHt.
$$

### Free Surfaces

The top and bottom edges are traction-free:

$$
\boldsymbol{\sigma}\mathbf{n}=\mathbf{0}
\qquad \text{at } y=0,H.
$$

## Finite-Element Discretization

Within each element:

$$
\mathbf{u}_h=\mathbf{N}\mathbf{d}_e,
\qquad
\boldsymbol{\varepsilon}_v=\mathbf{B}\mathbf{d}_e.
$$

Here, $\mathbf{d}_e$ contains element nodal displacements,
$\mathbf{N}$ contains shape functions, and $\mathbf{B}$
contains their spatial derivatives.

The element stiffness matrix is:

$$
\mathbf{K}_e
=
t\int_{\Omega_e}
\mathbf{B}^{T}\mathbf{D}\mathbf{B}\,dA.
$$

The consistent element load vector on a loaded edge is:

$$
\mathbf{F}_e
=
t\int_{\Gamma_{t,e}}
\mathbf{N}^{T}\overline{\mathbf{t}}\,ds.
$$

After assembly and application of the fixed-edge constraints:

$$
\mathbf{K}_{ff}\mathbf{d}_f=\mathbf{F}_f.
$$

Element stresses are recovered using:

$$
\boldsymbol{\sigma}_v
=
\mathbf{D}\mathbf{B}\mathbf{d}_e.
$$

Strain and stress may vary within an element, depending on
the element type.

## Visualization

The displayed deformed coordinates are:

$$
\mathbf{x}_{\mathrm{plot}}
=
\mathbf{x}+s\mathbf{u}_h.
$$

Here, $s$ is a visualization scale factor.
Use $s=1$ for the actual computed displacement.

The plane-stress von Mises equivalent stress is:

$$
\sigma_{\mathrm{vm}}
=
\sqrt{
\sigma_{xx}^{2}
-\sigma_{xx}\sigma_{yy}
+\sigma_{yy}^{2}
+3\tau_{xy}^{2}
}.
$$

## Verification

For a slender cantilever, Euler–Bernoulli beam theory gives
the approximate tip-deflection magnitude:

$$
\delta_{\mathrm{tip}}
\approx
\frac{PL^{3}}{3EI},
\qquad
I=\frac{tH^{3}}{12}.
$$

This is an approximate comparison for the 2D elasticity
solution; local end effects and shear deformation can
cause differences.

Check that support reactions balance the applied force.
Refine the mesh and monitor displacement convergence.

Peak stresses near the idealized fixed-edge corners may
be singular and should not be the sole convergence metric.

## Independent Variables

| File | Inputs |
| --- | --- |
| `geometry.m` | Length, height, thickness, and mesh definition |
| `material.m` | Young's modulus, Poisson's ratio, and constitutive assumption |
| `loads.m` | End traction or equivalent nodal forces |

These entries describe intended responsibilities; confirm
the actual variable names in the source files.

## How to Run

In MATLAB:

```matlab
main