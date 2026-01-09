function [nodes, elements] = generate_mesh_2D(Lx, Ly, nx, ny)
% Generates a structured rectangular mesh
dx = Lx/nx; dy = Ly/ny;
[x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);
nodes = [x(:), y(:)];  % Node coordinates

elements = zeros(nx*ny,4); % 4-node quad elements
count = 1;
for j = 1:ny
    for i = 1:nx
        n1 = (j-1)*(nx+1)+i;
        n2 = n1 + 1;
        n3 = n2 + nx +1;
        n4 = n1 + nx +1;
        elements(count,:) = [n1 n2 n3 n4];
        count = count +1;
    end
end
end
