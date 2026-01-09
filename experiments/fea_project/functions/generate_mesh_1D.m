function [nodes, elements] = generate_mesh_1D(L, n)
    nodes = linspace(0, L, n+1)';        % Node positions
    elements = [(1:n)' (2:n+1)'];        % Element connectivity
end
