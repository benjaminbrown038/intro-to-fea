# fea 3d project 

A MATLAB-based finite element analysis project that extends the FEA workflow to three-dimensional solid mechanics using tetrahedral elements. It solves for full 3D nodal displacements (uₓ, uᵧ, u_z), computes 3D strain and stress tensors, and visualizes results such as von Mises stress on volumetric meshes.

This project introduces realistic 3D geometries, material properties, and boundary conditions, forming the baseline for advanced capabilities like thermal–mechanical coupling, adaptive meshing, multi-material modeling, and industrial-style post-processing.

```
cd path_to/FEA_3D_Project
addpath(genpath(pwd))
main

```