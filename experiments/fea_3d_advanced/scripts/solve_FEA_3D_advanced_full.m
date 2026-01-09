run('problem_setup_3D_advanced.m');

[nodes, elements, element_material] = generate_mesh_3D(Lx,Ly,Lz,nx,ny,nz);

n_nodes = size(nodes,1);
Fvec = zeros(3*n_nodes,1);

% Apply mechanical load at top-front-right
load_node = find(nodes(:,1)==Lx & nodes(:,2)==Ly & nodes(:,3)==Lz);

fixed_nodes = find(nodes(:,3)==0);

for iter = 1:max_iterations
    fprintf('Adaptive iteration %d\n', iter);
    
    % Transient thermal loop
    for t = 1:nt
        DeltaT_elem = Tmax * sin(pi*t/nt) * mean(nodes(elements,3),2)/max(nodes(:,3));
        K = assemble_stiffness_3D(materials, element_material, nodes, elements);
        
        % Apply BCs
        Fmod = Fvec;
        for i=1:length(fixed_nodes)
            dof = 3*fixed_nodes(i)-2:3*fixed_nodes(i);
            K(dof,:) = 0; K(:,dof)=0; K(dof,dof)=eye(3); Fmod(dof)=0;
        end
        Fmod(3*load_node) = -F; % apply mechanical load

        % Solve linear or nonlinear iteration for plasticity
        U = K \ Fmod;
        
        [strain, stress] = compute_strain_stress_3D(U, nodes, elements, materials, element_material, DeltaT_elem);
        
        % Compute von Mises stress
        sigma_vm = sqrt(stress(1,:).^2 + stress(2,:).^2 + stress(3,:).^2 - ...
                        stress(1,:).*stress(2,:) - stress(2,:).*stress(3,:) - stress(3,:).*stress(1,:) + ...
                        3*(stress(4,:).^2 + stress(5,:).^2 + stress(6,:).^2));
        
        % Auto save figures per timestep
        plot_results_3D_advanced(nodes, elements, stress, U, scale_factor);
        saveas(gcf, fullfile(results_dir, sprintf('FEA_iter%d_t%d.png', iter, t)));
        
        % Export for Paraview
        export_vtk(nodes, elements, sigma_vm, U, fullfile(results_dir,sprintf('FEA_iter%d_t%d', iter,t)));
    end
    
    % Adaptive refinement (fully conforming)
    [nodes, elements, element_material] = adaptive_refine_3D_conforming(nodes, elements, element_material, sigma_vm, refine_threshold);
end

% Optional: export STL for CAD
export_stl(nodes, elements, fullfile(results_dir,'FEA_final_geometry'));
