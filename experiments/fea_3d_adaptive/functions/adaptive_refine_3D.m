function [new_nodes, new_elements] = adaptive_refine_3D(nodes, elements, sigma_vm, threshold)
max_sigma = max(sigma_vm);
high_stress_elems = find(sigma_vm > threshold*max_sigma);

new_nodes = nodes;
new_elements = elements;

for e = high_stress_elems
    elem_nodes = elements(e,:);
    % Create midpoints
    mid = zeros(6,3); idx = size(new_nodes,1)+1;
    edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for i=1:6
        mid(i,:) = mean(nodes(elem_nodes(edge_pairs(i,:)),:),1);
        new_nodes = [new_nodes; mid(i,:)]; %#ok<AGROW>
    end
    m = idx:idx+5; n = elem_nodes;
    % Replace tetrahedron with 8 smaller tetrahedra (conceptual)
    new_tets = [
        n(1) m(1) m(2) m(3);
        m(1) n(2) m(4) m(5);
        m(2) m(4) n(3) m(6);
        m(3) m(5) m(6) n(4);
        m(1) m(2) m(3) m(6);
        m(1) m(4) m(5) m(6);
        m(2) m(4) m(5) m(6);
        m(3) m(4) m(5) m(6);
    ];
    new_elements(e,:) = [];
    new_elements = [new_elements; new_tets];
end
end
