# fea project

A foundational MATLAB finite element analysis framework that implements the core FEA workflow: geometry definition, mesh generation, stiffness matrix assembly, boundary condition application, solving for nodal displacements, and post-processing of strain and stress.

It is designed as a learning and extensible base project, typically starting with simple linear elastic problems (often 1D or very simple 2D), and serves as the backbone from which the more advanced 2D, triangular, adaptive, and 3D FEA projects are built.

```
cd path_to/FEA_Project
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

