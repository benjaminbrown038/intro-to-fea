run('problem_setup_3D.m');

[nodes, elements] = generate_mesh_3D(Lx,Ly,Lz,nx,ny,nz);

K = assemble_stiffness_3D(E, nu, nodes, elements);

n_nodes = size(nodes,1);
Fvec = zeros(3*n_nodes,1);

% Apply mechanical load at top right front node
load_node = find(nodes(:,1)==Lx & nodes(:,2)==Ly & nodes(:,3)==Lz);
Fvec(3*load_node) = -F;  % z-direction

% Fix bottom layer nodes (z=0)
fixed_nodes = find(nodes(:,3)==0);
for i=1:length(fixed_nodes)
    dof = 3*fixed_nodes(i)-2:3*fixed_nodes(i);
    K(dof,:) = 0; K(:,dof)=0; K(dof,dof)=eye(3); Fvec(dof)=0;
end

U = K \ Fvec;

[strain, stress] = compute_strain_stress_3D(U, nodes, elements, E, nu, alpha, Tmax);

% Simple visualization (slice or scatter)
sigma_vm = sqrt(stress(1,:).^2 + stress(2,:).^2 + stress(3,:).^2 - stress(1,:).*stress(2,:) - stress(2,:).*stress(3,:) - stress(3,:).*stress(1,:) + 3*(stress(4,:).^2 + stress(5,:).^2 + stress(6,:).^2));
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 20, sigma_vm, 'filled'); colorbar;
title('3D Von Mises Stress'); axis equal; view(3);
