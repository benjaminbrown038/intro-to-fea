function [nodes, elements] = generate_tri_mesh(geom)

[x,y] = meshgrid(0:geom.h:geom.L, 0:geom.h:geom.H);
nodes = [x(:), y(:)];

DT = delaunayTriangulation(nodes);
elements = DT.ConnectivityList;
end
