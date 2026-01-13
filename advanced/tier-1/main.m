clear; clc; close all;

% Load parameters
geom = geometry();
mat  = material();
load = loads();

% Generate mesh
[nodes, elements] = generate_mesh(geom);

% Assemble global stiffness matrix
K = assemble_global_K(nodes, elements, mat);

% Apply boundary conditions and loads
[F, fixed_dofs] = apply_boundary_conditions(nodes, load, geom);

% Solve system
U = solve_system(K, F, fixed_dofs);

% Post-processing
[strain, stress] = compute_strain_stress(nodes, elements, U, mat);
plot_results(nodes, elements, U, stress);
