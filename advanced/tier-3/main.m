clear; clc; close all;

geom = geometry();
mat  = material();
load = loads();

[nodes, elements] = generate_tri_mesh(geom);

K = assemble_global_K(nodes, elements, mat);

[F, fixed_dofs] = apply_boundary_conditions(nodes, load, geom);

U = solve_system(K, F, fixed_dofs);

[strain, stress] = compute_strain_stress(nodes, elements, U, mat);

plot_results(nodes, elements, U, stress);
