function [strain, stress] = compute_strain_stress(nodes, elements, U, mat)

ne = size(elements,1);
strain = zeros(ne,3);
stress = zeros(ne,3);

for e = 1:ne
    elem = elements(e,:);
    coords = nodes(elem,:);
    dofs = reshape([2*elem-1; 2*elem],1,[]);
    ue = U(dofs);

    x = coords(:,1); y = coords(:,2);
    A = polyarea(x,y);

    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];

    B = 1/(2*A) * ...
        [b(1) 0 b(2) 0 b(3) 0;
         0 c(1) 0 c(2) 0 c(3);
         c(1) b(1) c(2) b(2) c(3) b(3)];

    strain(e,:) = B * ue;
    stress(e,:) = mat.D * strain(e,:)';
end
end
