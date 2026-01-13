function [strain, stress] = compute_strain_stress(nodes, elements, U, mat)

D = (mat.E/(1-mat.nu^2))*[1 mat.nu 0; mat.nu 1 0; 0 0 (1-mat.nu)/2];

stress = zeros(size(elements,1),3);
strain = zeros(size(elements,1),3);

for e = 1:size(elements,1)
    elem_nodes = elements(e,:);
    coords = nodes(elem_nodes,:);
    dofs = reshape([2*elem_nodes-1; 2*elem_nodes],1,[]);
    ue = U(dofs);

    B = zeros(3,8);
    B(:,1:2:end) = [1 0; 0 0; 0 1];
    B(:,2:2:end) = [0 0; 0 1; 1 0];

    strain(e,:) = B*ue;
    stress(e,:) = D*strain(e,:)';
end
end
