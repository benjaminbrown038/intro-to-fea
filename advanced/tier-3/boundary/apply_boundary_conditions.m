function [F, fixed_dofs] = apply_boundary_conditions(nodes, load, geom)

ndof = size(nodes,1)*2;
F = zeros(ndof,1);

right = find(abs(nodes(:,1)-geom.L) < 1e-6);
F(2*right) = load.Fy / length(right);

left = find(abs(nodes(:,1)) < 1e-6);
fixed_dofs = reshape([2*left-1; 2*left],1,[]);
end
