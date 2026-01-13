function K = assemble_global_K(nodes, elements, mat)

ndof = size(nodes,1)*2;
K = zeros(ndof);

for e = 1:size(elements,1)
    elem_nodes = elements(e,:);
    coords = nodes(elem_nodes,:);
    Ke = element_stiffness(mat.E, mat.nu, 1.0, coords);

    dofs = reshape([2*elem_nodes-1; 2*elem_nodes],1,[]);
    K(dofs,dofs) = K(dofs,dofs) + Ke;
end
end
