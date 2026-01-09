function plot_results_3D(nodes, elements, sigma_vm, U, scale_factor)
% 3D FEA results visualization
% nodes     : [n x 3] coordinates
% elements  : [m x 4] tetrahedral connectivity
% sigma_vm  : [1 x m] von Mises stress per element
% U         : [3*n x 1] nodal displacement
% scale_factor : exaggeration factor for deformation

if nargin < 5
    scale_factor = 1; % default no exaggeration
end

% Deformed nodal coordinates
nodes_def = nodes + scale_factor * reshape(U,3,[]);

% Plot color-mapped tetrahedra
figure('Color','w'); hold on;

for e = 1:size(elements,1)
    tet_nodes = elements(e,:);
    coords = nodes_def(tet_nodes,:);
    patch('Vertices', coords, 'Faces', [1 2 3 4],
          'FaceColor','interp','FaceVertexCData', sigma_vm(e)*ones(4,1),
          'EdgeColor','k','LineWidth',0.5);
end

colormap(jet); colorbar; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Adaptive FEA: Von Mises Stress with Deformation');

% Optional: slice plane along Z = mid
z_slice = mean(nodes(:,3));
idx = find(abs(nodes(:,3) - z_slice) < 1e-6);
scatter3(nodes_def(idx,1), nodes_def(idx,2), nodes_def(idx,3), 50, sigma_vm(idx), 'filled');
end
