# fea 2d project

A MATLAB-based 2D finite element analysis solver for planar solid mechanics.
It models structures using 2D elements (triangular or quadrilateral) under plane stress or plane strain assumptions, solves for nodal displacements (uₓ, uᵧ), computes element strains and stresses, and visualizes results using stress contour plots (e.g., von Mises).

Typical uses: plates, brackets, cross-sections, and validation of FEA fundamentals before moving to 3D.


```
cd path_to/FEA_2D_Project
addpath(genpath(pwd))
main

```

### Geometry 

Length 
Height
Thickness
Cross Section

### Material Properties 

Young's Modulus 
Poisson's Ratio
Density 
Thermal Expansion


### Boundary Conditions

Fixed Edge
Roller Edge
Symmetry
Custom Constraints


### Mesh Parameters

Elements Along X
Elements Along Y 
Element Type 
Mesh Density 


### Loading 

Load Type 
Load Magnitude 
Load Direction
Time Dependence

### Simulation Options 

Solver Type 
Analysis Type 
Linear/Nonlinear 
Post Processing

