function [nodes, elements] = generate_mesh_3D(Lx,Ly,Lz,nx,ny,nz)
% Generate a simple 3D grid and tetrahedral mesh using delaunay

[x, y, z] = ndgrid(linspace(0,Lx,nx+1), linspace(0,Ly,ny+1), linspace(0,Lz,nz+1));
nodes = [x(:), y(:), z(:)];

DT = delaunayTriangulation(nodes);
elements = DT.ConnectivityList; % 4-node tetrahedra
end
