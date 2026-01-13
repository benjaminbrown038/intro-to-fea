function [nodes, elements] = generate_mesh(geom)

[x,y] = meshgrid( ...
    linspace(0,geom.L,geom.nx+1), ...
    linspace(0,geom.H,geom.ny+1));

nodes = [x(:), y(:)];

elements = [];
nx = geom.nx;

for j = 1:geom.ny
    for i = 1:geom.nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + nx + 1;
        n4 = n3 - 1;
        elements = [elements; n1 n2 n3 n4];
    end
end
end
