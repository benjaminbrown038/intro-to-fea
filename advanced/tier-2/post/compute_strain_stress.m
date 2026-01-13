function [strain, stress] = compute_strain_stress(nodes, elements, U, mat)

ne = size(elements,1);
strain = zeros(ne,3);
stress = zeros(ne,3);

for e = 1:ne
    elem = elements(e,:);
    coords = nodes(elem,:);
    dofs = reshape([2*elem-1; 2*elem],1,[]);
    ue = U(dofs);

    % Evaluate at element center
    xi = 0; eta = 0;

    dNdxi = 0.25 * [-(1-eta) (1-eta) (1+eta) -(1+eta)];
    dNdeta= 0.25 * [-(1-xi) -(1+xi) (1+xi) (1-xi)];

    J = [dNdxi; dNdeta]*coords;
    dN = inv(J)*[dNdxi; dNdeta];

    B = zeros(3,8);
    for k = 1:4
        B(:,2*k-1:2*k) = ...
            [dN(1,k) 0;
             0 dN(2,k);
             dN(2,k) dN(1,k)];
    end

    strain(e,:) = B*ue;
    stress(e,:) = mat.D * strain(e,:)';
end
end
