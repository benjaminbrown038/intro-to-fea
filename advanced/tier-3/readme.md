# Tier 3 – 2D Triangular Element FEA

## Overview

Analyze 2D linear elastic structures using three-node
constant-strain triangular elements (CST).

Triangular meshes accommodate irregular boundaries and
nonuniform element sizes. Curved boundaries are approximated
by straight element edges.

## Physical Assumptions

- Static equilibrium.
- Small strains and small rotations.
- Linear elastic, isotropic material.
- Plane stress.
- Constant out-of-plane thickness.
- Material properties constant within each element.

## Governing Equations

### Static Equilibrium

$$
-\nabla\cdot\boldsymbol{\sigma}=\mathbf{b}
\qquad \text{in } \Omega.
$$

Here, $\mathbf{b}$ is body force per unit volume.

For problems loaded only through boundary tractions:

$$
\mathbf{b}=\mathbf{0}.
$$

### Strain–Displacement Relationship

Using engineering shear strain:

$$
\boldsymbol{\varepsilon}_v=
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

Engineering shear strain satisfies:

$$
\gamma_{xy}=2\varepsilon_{xy}.
$$

### Plane-Stress Material Model

$$
\boldsymbol{\sigma}_v
=
\mathbf{D}\boldsymbol{\varepsilon}_v,
\qquad
\boldsymbol{\sigma}_v=
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

## Triangular Element Geometry

Each triangle contains three nodes:

$$
(x_1,y_1),\qquad
(x_2,y_2),\qquad
(x_3,y_3).
$$

Order the nodes counterclockwise. The element area is:

$$
A_e=
\frac{1}{2}
\det
\begin{bmatrix}
1 & x_1 & y_1\\
1 & x_2 & y_2\\
1 & x_3 & y_3
\end{bmatrix}
>0.
$$

Reorder clockwise elements before evaluating the formulas.
Reject zero-area or nearly degenerate triangles.

## Linear Shape Functions

For cyclic node indices
$(i,j,k)=(1,2,3),(2,3,1),(3,1,2)$, define:

$$
a_i=x_jy_k-x_ky_j,
$$

$$
b_i=y_j-y_k,
\qquad
c_i=x_k-x_j.
$$

The shape functions are:

$$
N_i(x,y)=\frac{a_i+b_ix+c_iy}{2A_e}.
$$

They satisfy:

$$
N_i(x_j,y_j)=\delta_{ij},
\qquad
N_1+N_2+N_3=1.
$$

## Displacement Interpolation

Each node has two displacement degrees of freedom:

$$
\mathbf{d}_e=
\begin{bmatrix}
u_{x,1}\\
u_{y,1}\\
u_{x,2}\\
u_{y,2}\\
u_{x,3}\\
u_{y,3}
\end{bmatrix}.
$$

The displacement inside the triangle is:

$$
\mathbf{u}_h=\mathbf{N}\mathbf{d}_e,
$$

where:

$$
\mathbf{N}=
\begin{bmatrix}
N_1 & 0 & N_2 & 0 & N_3 & 0\\
0 & N_1 & 0 & N_2 & 0 & N_3
\end{bmatrix}.
$$

## Constant Strain–Displacement Matrix

Because the shape functions are linear:

$$
\frac{\partial N_i}{\partial x}
=
\frac{b_i}{2A_e},
\qquad
\frac{\partial N_i}{\partial y}
=
\frac{c_i}{2A_e}.
$$

Therefore:

$$
\mathbf{B}
=
\frac{1}{2A_e}
\begin{bmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0\\
0 & c_1 & 0 & c_2 & 0 & c_3\\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{bmatrix}.
$$

The element strain is:

$$
\boldsymbol{\varepsilon}_v
=
\mathbf{B}\mathbf{d}_e.
$$

Since $\mathbf{B}$ is constant, strain is constant throughout
each triangle.

## Element Stiffness Matrix

For constant thickness $t$ and constitutive matrix $\mathbf{D}$:

$$
\mathbf{K}_e
=
t\int_{\Omega_e}
\mathbf{B}^{T}\mathbf{D}\mathbf{B}\,dA.
$$

This reduces exactly to:

$$
\mathbf{K}_e
=
tA_e\mathbf{B}^{T}\mathbf{D}\mathbf{B}.
$$

Each element contributes a $6\times6$ stiffness matrix.

## Element Loads

### Constant Body Force

For a constant body-force vector
$\mathbf{b}=[b_x,b_y]^T$:

$$
\mathbf{F}_e^{\mathrm{body}}
=
\frac{tA_e}{3}
\begin{bmatrix}
b_x\\
b_y\\
b_x\\
b_y\\
b_x\\
b_y
\end{bmatrix}.
$$

### Uniform Edge Traction

For a straight edge of length $\ell$ connecting nodes $i$
and $j$, a constant traction $\overline{\mathbf{t}}$
produces equal consistent nodal loads:

$$
\mathbf{F}_i^{\mathrm{edge}}
=
\mathbf{F}_j^{\mathrm{edge}}
=
\frac{t\ell}{2}\overline{\mathbf{t}}.
$$

Traction is force per unit area; the factor $t\ell$
converts it to force.

## Boundary Conditions

Prescribed displacement:

$$
\mathbf{u}=\overline{\mathbf{u}}
\qquad \text{on } \Gamma_D.
$$

Prescribed traction:

$$
\boldsymbol{\sigma}\mathbf{n}
=
\overline{\mathbf{t}}
\qquad \text{on } \Gamma_N.
$$

A fixed boundary has $\overline{\mathbf{u}}=\mathbf{0}$.
A free boundary has $\overline{\mathbf{t}}=\mathbf{0}$.

## Global Assembly and Solution

Let $\mathbf{A}_e$ extract element degrees of freedom:

$$
\mathbf{d}_e=\mathbf{A}_e\mathbf{d}.
$$

Assemble:

$$
\mathbf{K}
=
\sum_e\mathbf{A}_e^T\mathbf{K}_e\mathbf{A}_e,
$$

$$
\mathbf{F}
=
\sum_e\mathbf{A}_e^T\mathbf{F}_e.
$$

Separate free and prescribed degrees of freedom:

$$
\mathbf{K}_{ff}\mathbf{d}_f
=
\mathbf{F}_f-\mathbf{K}_{fc}\mathbf{d}_c.
$$

Supports must remove unconstrained rigid-body motion.

## Stress Recovery

Element stress is:

$$
\boldsymbol{\sigma}_{v,e}
=
\mathbf{D}\mathbf{B}\mathbf{d}_e.
$$

For a homogeneous CST element, stress is constant within
the triangle but may differ between neighboring triangles.

Plane-stress von Mises stress is:

$$
\sigma_{\mathrm{vm},e}
=
\sqrt{
\sigma_{xx,e}^{2}
-\sigma_{xx,e}\sigma_{yy,e}
+\sigma_{yy,e}^{2}
+3\tau_{xy,e}^{2}
}.
$$

## Visualization

The displayed deformed coordinates are:

$$
\mathbf{x}_{\mathrm{plot}}
=
\mathbf{x}+s\mathbf{u}_h.
$$

Here, $s$ is a display scale factor.

Use one color per triangle to display the computed
element stress directly.

If nodal averaging or interpolation is used to produce
smooth contours, label the result as recovered or smoothed
stress. Avoid averaging across different materials.

## Mesh Generation

Delaunay triangulation can generate triangle connectivity
from a set of points.

For nonconvex boundaries or holes, use boundary constraints
and domain filtering as appropriate. Unconstrained Delaunay
triangulation alone does not identify the intended domain.

Check:

- Positive element areas.
- Consistent node ordering.
- Correct representation of boundaries and holes.
- No overlapping or degenerate elements.
- Adequate resolution near high stress gradients.

Triangular elements support nonuniform meshes, but adaptive
refinement also requires an error indicator, element marking,
and a conforming refinement procedure.

## Verification

### Constant-Strain Patch Test

Apply an affine displacement field:

$$
u_x=\alpha_0+\alpha_1x+\alpha_2y,
$$

$$
u_y=\beta_0+\beta_1x+\beta_2y.
$$

The exact strain is constant:

$$
\boldsymbol{\varepsilon}_v=
\begin{bmatrix}
\alpha_1\\
\beta_2\\
\alpha_2+\beta_1
\end{bmatrix}.
$$

A correctly implemented CST mesh should reproduce this
field under compatible prescribed boundary displacements,
up to numerical roundoff.

### Additional Checks

- Rigid translation produces zero strain.
- Infinitesimal rigid rotation produces zero strain.
- Support reactions balance external loads.
- Displacements converge under mesh refinement.

CST elements can be overly stiff in bending on coarse
meshes. Use sufficient refinement, especially through
the thickness of a bending region.

## How to Run

In MATLAB:

```matlab
main
```

## Outputs

- Two displacement components at each node.
- Constant strain per triangle.
- Constant stress per triangle.
- Element stress contour plots.
- Deformed mesh visualization.