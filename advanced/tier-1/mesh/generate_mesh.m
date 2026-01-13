function [nodes, elements] = generate_mesh(geom)

nx = geom.nx;
ny = geom.ny;
L  = geom.L;
H  = geom.H;

[x, y] = meshgrid(linspace(0,L,nx+1), linspace(0,H,ny+1));
nodes = [x(:), y(:)];

elements = [];
for j = 1:ny
    for i = 1:nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + nx + 1;
        n4 = n3 - 1;
        elements = [elements; n1 n2 n3 n4];
    end
end
end
