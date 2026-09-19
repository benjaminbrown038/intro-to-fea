# Tier 2 – 2D Linear Elastic FEA

## Overview

Analyze a plane-stress cantilever beam using a structured
mesh of four-node rectangular elements.

This tier emphasizes element stiffness matrices, global
assembly, 2D displacement fields, and stress visualization.

## Physical Assumptions

- Static equilibrium.
- Small strains and small rotations.
- Homogeneous, isotropic, linear elastic material.
- Plane stress with constant out-of-plane thickness.
- Fixed left edge and prescribed loading on the right edge.

## Geometry and Mesh

The domain is:

$$
\Omega=[0,L]\times[0,H].
$$

For a uniform mesh with $n_x$ elements along the length
and $n_y$ along the height:

$$
h_x=\frac{L}{n_x},
\qquad
h_y=\frac{H}{n_y}.
$$

For four-node rectangular elements:

$$
N_{\mathrm{elements}}=n_xn_y,
$$

$$
N_{\mathrm{nodes}}=(n_x+1)(n_y+1).
$$

Each node has two displacement degrees of freedom:

$$
\mathbf{d}_i=
\begin{bmatrix}
u_{x,i}\\
u_{y,i}
\end{bmatrix}.
$$

Before applying constraints:

$$
N_{\mathrm{DOF}}=2N_{\mathrm{nodes}}.
$$

## Governing Equations

### Static Equilibrium

$$
-\nabla\cdot\boldsymbol{\sigma}=\mathbf{b}.
$$

Here, $\mathbf{b}$ is body force per unit volume.
For loading applied only at the beam end:

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
\begin{bmatrix}
\frac{\partial u_x}{\partial x}\\
\frac{\partial u_y}{\partial y}\\
\frac{\partial u_x}{\partial y}
\frac{\partial u_y}{\partial x}
\end{bmatrix}.
$$

### Plane-Stress Material Model

$$
\boldsymbol{\sigma}_v
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

## Four-Node Quadrilateral Element

### Natural Coordinates

The reference element uses:

$$
-1\le\xi\le1,
\qquad
-1\le\eta\le1.
$$

Nodes are numbered counterclockwise:

$$
(\xi_i,\eta_i)
(-1,-1),\ (1,-1),\ (1,1),\ (-1,1).
$$

### Shape Functions

$$
N_1=\frac{1}{4}(1-\xi)(1-\eta),
\qquad
N_2=\frac{1}{4}(1+\xi)(1-\eta),
$$

$$
N_3=\frac{1}{4}(1+\xi)(1+\eta),
\qquad
N_4=\frac{1}{4}(1-\xi)(1+\eta).
$$

Geometry and displacement use the same interpolation:

$$
x=\sum_{i=1}^{4}N_i x_i,
\qquad
y=\sum_{i=1}^{4}N_i y_i,
$$

$$
u_x=\sum_{i=1}^{4}N_i u_{x,i},
\qquad
u_y=\sum_{i=1}^{4}N_i u_{y,i}.
$$

### Coordinate Transformation

Define the Jacobian:

$$
\mathbf{J}
\begin{bmatrix}
\frac{\partial x}{\partial\xi}
\frac{\partial x}{\partial\eta}\\
\frac{\partial y}{\partial\xi}
\frac{\partial y}{\partial\eta}
\end{bmatrix}.
$$

Shape-function derivatives in physical coordinates are:

$$
\begin{bmatrix}
\frac{\partial N_i}{\partial x}\\
\frac{\partial N_i}{\partial y}
\end{bmatrix}
\mathbf{J}^{-T}
\begin{bmatrix}
\frac{\partial N_i}{\partial\xi}\\
\frac{\partial N_i}{\partial\eta}
\end{bmatrix}.
$$

For a valid, consistently oriented element:

$$
\det(\mathbf{J})>0.
$$

### Strain–Displacement Matrix

For node $i$:

$$
\mathbf{B}_i=
\begin{bmatrix}
\frac{\partial N_i}{\partial x} & 0\\
0 & \frac{\partial N_i}{\partial y}\\
\frac{\partial N_i}{\partial y}
&
\frac{\partial N_i}{\partial x}
\end{bmatrix}.
$$

The complete element matrix is:

$$
\mathbf{B}
\begin{bmatrix}
\mathbf{B}_1 & \mathbf{B}_2 &
\mathbf{B}_3 & \mathbf{B}_4
\end{bmatrix}.
$$

With nodal ordering:

$$
\mathbf{d}_e=
\begin{bmatrix}
u_{x,1}&u_{y,1}&
u_{x,2}&u_{y,2}&
u_{x,3}&u_{y,3}&
u_{x,4}&u_{y,4}
\end{bmatrix}^{T},
$$

the element strain is:

$$
\boldsymbol{\varepsilon}_v=\mathbf{B}\mathbf{d}_e.
$$

## Element Stiffness

For thickness $t$:

$$
\mathbf{K}_e
t\int_{-1}^{1}\int_{-1}^{1}
\mathbf{B}^{T}\mathbf{D}\mathbf{B}
\det(\mathbf{J})
\,d\xi\,d\eta.
$$

Using standard $2\times2$ Gauss integration:

$$
\mathbf{K}_e
\approx
t\sum_{g=1}^{4}
\mathbf{B}_g^{T}\mathbf{D}\mathbf{B}_g
\det(\mathbf{J}_g)w_g.
$$

The integration points are all combinations of:

$$
\xi_g,\eta_g=\pm\frac{1}{\sqrt{3}},
\qquad
w_g=1.
$$

## Loads and Boundary Conditions

### Fixed Edge

At $x=0$:

$$
u_x=u_y=0.
$$

### Loaded Edge

At $x=L$:

$$
\boldsymbol{\sigma}\mathbf{n}
\overline{\mathbf{t}}.
$$

For uniform downward traction of magnitude $q$:

$$
\overline{\mathbf{t}}
\begin{bmatrix}
0\\
-q
\end{bmatrix},
\qquad
P=qHt.
$$

Here, $P$ is the total downward force magnitude.

The consistent edge-load vector is:

$$
\mathbf{F}_e^{\mathrm{edge}}
t\int_{\Gamma_{t,e}}
\mathbf{N}^{T}\overline{\mathbf{t}}\,ds.
$$

The top and bottom edges are traction-free.

## Global Assembly

Let $\mathbf{A}_e$ extract an element's degrees of freedom
from the global displacement vector:

$$
\mathbf{d}_e=\mathbf{A}_e\mathbf{d}.
$$

Assemble the global stiffness and load:

$$
\mathbf{K}
\sum_e\mathbf{A}_e^{T}\mathbf{K}_e\mathbf{A}_e,
$$

$$
\mathbf{F}
\sum_e\mathbf{A}_e^{T}\mathbf{F}_e.
$$

After applying zero displacement at the fixed edge:

$$
\mathbf{K}_{ff}\mathbf{d}_f=\mathbf{F}_f.
$$

Using the original assembled matrices, support reactions are:

$$
\mathbf{R}_c
\mathbf{K}_{cf}\mathbf{d}_f-\mathbf{F}_c.
$$

## Strain and Stress Recovery

At each evaluation point:

$$
\boldsymbol{\varepsilon}_v
\mathbf{B}\mathbf{d}_e,
$$

$$
\boldsymbol{\sigma}_v
\mathbf{D}\mathbf{B}\mathbf{d}_e.
$$

For Q4 elements, strain and stress generally vary within
each element.

Plane-stress von Mises stress is:

$$
\sigma_{\mathrm{vm}}
\sqrt{
\sigma_{xx}^{2}
-\sigma_{xx}\sigma_{yy}
+\sigma_{yy}^{2}
+3\tau_{xy}^{2}
}.
$$

State whether contour values are evaluated at element
centers, integration points, or recovered nodes.
Nodal averaging can smooth discontinuities between elements.

## Deformed Shape

$$
\mathbf{x}_{\mathrm{plot}}
\mathbf{x}+s\mathbf{u}_h.
$$

Here, $s$ is the display scale factor; $s=1$ shows the actual
computed displacement.

## Verification

For a slender beam, compare the tip-deflection magnitude
with the Euler–Bernoulli approximation:

$$
\delta_{\mathrm{tip}}
\approx
\frac{PL^3}{3EI},
\qquad
I=\frac{tH^3}{12}.
$$

Also check:

- Support reactions balance the applied force and moment.
- Displacement changes decrease as the mesh is refined.
- Element Jacobian determinants remain positive.
- Stress contours identify their evaluation or recovery method.

Fully integrated Q4 elements can be excessively stiff in
bending on coarse meshes. Refine through the beam height
as well as along its length.

Idealized fixed-edge corner stresses may be singular;
use displacement and stresses away from those corners
when assessing convergence.

## Independent Variables

| Input | Location |
| --- | --- |
| Geometry and thickness | `geometry.m` |
| Young's modulus and Poisson's ratio | `material.m` |
| End loading | `loads.m` |
| Mesh density | `nx`, `ny` in the mesh setup |

The mesh formulas above treat `nx` and `ny` as element
counts; confirm that convention in the source code.

## How to Run

In MATLAB:

```matlab
main
```

## Outputs

- Two displacement components at each node.
- Element strain and stress.
- Stress contour plots.
- Deformed shape.