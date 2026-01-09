run('problem_setup_3D_adaptive.m');

[nodes, elements] = generate_mesh_3D(Lx,Ly,Lz,nx,ny,nz);

% Force vector
n_nodes = size(nodes,1);
Fvec = zeros(3*n_nodes,1);
load_node = find(nodes(:,1)==Lx & nodes(:,2)==Ly & nodes(:,3)==Lz);
Fvec(3*load_node) = -F;

% Fix bottom layer
fixed_nodes = find(nodes(:,3)==0);
for i=1:length(fixed_nodes)
    dof = 3*fixed_nodes(i)-2:3*fixed_nodes(i);
    K(dof,:) = 0; K(:,dof)=0; K(dof,dof)=eye(3); Fvec(dof)=0;
end

% Adaptive refinement loop
for iter = 1:max_iterations
    K = assemble_stiffness_3D(E, nu, nodes, elements);
    U = K \ Fvec;
    [strain, stress] = compute_strain_stress_3D(U, nodes, elements, E, nu, alpha, Tmax);

    sigma_vm = sqrt(stress(1,:).^2 + stress(2,:).^2 + stress(3,:).^2 - ...
                    stress(1,:).*stress(2,:) - stress(2,:).*stress(3,:) - stress(3,:).*stress(1,:) + ...
                    3*(stress(4,:).^2 + stress(5,:).^2 + stress(6,:).^2));

    [nodes, elements] = adaptive_refine_3D(nodes, elements, sigma_vm, refine_threshold);
end

% Simple visualization
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 20, sigma_vm, 'filled');
colorbar; axis equal; view(3); title('3D Von Mises Stress (Adaptive)');
