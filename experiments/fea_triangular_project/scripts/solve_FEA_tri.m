%% Load Problem Setup
run('problem_setup_tri.m');

%% Generate Triangular Mesh
[nodes, elements] = generate_mesh_tri(Lx, Ly, nx, ny);

%% Assemble Global Stiffness
K = assemble_stiffness_tri(E, nu, nodes, elements);

%% Force Vector
n_nodes = size(nodes,1);
Fvec = zeros(2*n_nodes,1);

% Example: apply vertical force at top right node
top_right_node = find(nodes(:,1)==Lx & nodes(:,2)==Ly);
Fvec(2*top_right_node) = -F;

%% Apply Boundary Conditions (fix bottom edge)
fixed_nodes = find(nodes(:,2)==0);
for i = 1:length(fixed_nodes)
    dof_x = 2*fixed_nodes(i)-1; dof_y = 2*fixed_nodes(i);
    K(dof_x,:) = 0; K(:,dof_x)=0; K(dof_x,dof_x)=1; Fvec(dof_x)=0;
    K(dof_y,:) = 0; K(:,dof_y)=0; K(dof_y,dof_y)=1; Fvec(dof_y)=0;
end

%% Solve
U = K \ Fvec;

%% Compute Strain and Stress
[strain, stress] = compute_strain_stress_tri(U, nodes, elements, E, nu);

%% Plot
plot_results_tri(nodes, elements, stress);
