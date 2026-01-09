function plot_results_3D_advanced(nodes, elements, stress, U, scale_factor)
% Advanced 3D FEA visualization
% Shows deformed mesh + von Mises stress + slice planes + principal stress vectors
%
% nodes    : [n x 3] nodal coordinates
% elements : [m x 4] tetrahedral connectivity
% stress   : [6 x m] element stress (Voigt)
% U        : [3*n x 1] nodal displacements
% scale_factor : exaggeration factor

if nargin < 5
    scale_factor = 1;
end

nodes_def = nodes + scale_factor * reshape(U,3,[]);

% Compute von Mises stress
sigma_vm = sqrt(stress(1,:).^2 + stress(2,:).^2 + stress(3,:).^2 - ...
                stress(1,:).*stress(2,:) - stress(2,:).*stress(3,:) - stress(3,:).*stress(1,:) + ...
                3*(stress(4,:).^2 + stress(5,:).^2 + stress(6,:).^2));

figure('Color','w'); hold on;

% Plot tetrahedral mesh colored by von Mises stress
for e = 1:size(elements,1)
    tet_nodes = elements(e,:);
    coords = nodes_def(tet_nodes,:);
    patch('Vertices', coords, 'Faces', [1 2 3 4], ...
          'FaceColor','interp','FaceVertexCData', sigma_vm(e)*ones(4,1), ...
          'EdgeColor','k','LineWidth',0.5);
end
colormap(jet); colorbar; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Adaptive FEA: Von Mises Stress with Deformation');

% Slice planes at mid X, Y, Z
x_mid = mean(nodes(:,1)); y_mid = mean(nodes(:,2)); z_mid = mean(nodes(:,3));

scatter3(nodes_def(abs(nodes(:,1)-x_mid)<1e-3,1), ...
         nodes_def(abs(nodes(:,1)-x_mid)<1e-3,2), ...
         nodes_def(abs(nodes(:,1)-x_mid)<1e-3,3), 50, ...
         sigma_vm(abs(nodes(:,1)-x_mid)<1e-3), 'filled');

scatter3(nodes_def(abs(nodes(:,2)-y_mid)<1e-3,1), ...
         nodes_def(abs(nodes(:,2)-y_mid)<1e-3,2), ...
         nodes_def(abs(nodes(:,2)-y_mid)<1e-3,3), 50, ...
         sigma_vm(abs(nodes(:,2)-y_mid)<1e-3), 'filled');

scatter3(nodes_def(abs(nodes(:,3)-z_mid)<1e-3,1), ...
         nodes_def(abs(nodes(:,3)-z_mid)<1e-3,2), ...
         nodes_def(abs(nodes(:,3)-z_mid)<1e-3,3), 50, ...
         sigma_vm(abs(nodes(:,3)-z_mid)<1e-3), 'filled');

% Principal stress vectors (simplified)
for e = 1:size(elements,1)
    % Convert Voigt stress to tensor
    S = [stress(1,e) stress(4,e) stress(6,e);
         stress(4,e) stress(2,e) stress(5,e);
         stress(6,e) stress(5,e) stress(3,e)];
    [V,D] = eig(S); % principal directions & magnitudes
    center = mean(nodes_def(elements(e,:),:),1);
    quiver3(center(1), center(2), center(3), V(1,1)*0.01, V(2,1)*0.01, V(3,1)*0.01, 'r', 'LineWidth',1.2);
end

rotate3d on;
end
