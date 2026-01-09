function K = assemble_stiffness_1D(E, A, nodes, elements)
    n_nodes = length(nodes);
    K = zeros(n_nodes);
    
    for e = 1:size(elements,1)
        n1 = elements(e,1);
        n2 = elements(e,2);
        Le = nodes(n2) - nodes(n1);
        ke = (E*A/Le) * [1 -1; -1 1];
        K(n1:n2, n1:n2) = K(n1:n2, n1:n2) + ke;
    end
end
