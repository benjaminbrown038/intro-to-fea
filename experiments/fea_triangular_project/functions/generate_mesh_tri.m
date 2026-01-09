function [nodes, elements] = generate_mesh_tri(Lx, Ly, nx, ny)
% Generates a triangular mesh using Delaunay triangulation

% Structured grid points
[x, y] = meshgrid(linspace(0,Lx,nx+1), linspace(0,Ly,ny+1));
nodes = [x(:), y(:)];

% Delaunay triangulation
DT = delaunayTriangulation(nodes);
elements = DT.ConnectivityList;
end
