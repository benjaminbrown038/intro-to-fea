# fea 3d adaptive

A 3D MATLAB finite element analysis project that builds on the basic 3D solver by introducing adaptive mesh refinement. The mesh is automatically refined in regions with high stress, strain, or temperature gradients, improving accuracy while keeping computational cost manageable.

The project solves 3D linear elastic problems with tetrahedral elements, evaluates error or gradient indicators, refines the mesh locally, and re-solves iteratively. It mirrors industrial FEA practice and forms the basis for high-fidelity simulations involving thermal stress, complex geometries, and localized failure analysis.

```
cd path_to/FEA_3D_Adaptive
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

