function K = assemble_global_K(nodes, elements, mat)

ndof = size(nodes,1)*2;
K = zeros(ndof);

for e = 1:size(elements,1)
    elem = elements(e,:);
    coords = nodes(elem,:);
    Ke = element_Q4(coords, mat, 1.0);

    dofs = reshape([2*elem-1; 2*elem],1,[]);
    K(dofs,dofs) = K(dofs,dofs) + Ke;
end
end
