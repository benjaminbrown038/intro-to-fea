function [F, fixed_dofs] = apply_boundary_conditions(nodes, load, geom)

ndof = size(nodes,1)*2;
F = zeros(ndof,1);

% Apply force at right edge
right_nodes = find(abs(nodes(:,1)-geom.L) < 1e-6);
F(2*right_nodes) = load.Fy/length(right_nodes);

% Fix left edge
left_nodes = find(abs(nodes(:,1)) < 1e-6);
fixed_dofs = reshape([2*left_nodes-1; 2*left_nodes],1,[]);
end
