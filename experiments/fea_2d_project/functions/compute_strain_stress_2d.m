function [strain, stress] = compute_strain_stress_2D(U, nodes, elements, E, nu)
C = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
n_elem = size(elements,1);
strain = zeros(3,n_elem);
stress = zeros(3,n_elem);

for e = 1:n_elem
    elem_nodes = elements(e,:);
    xe = nodes(elem_nodes,1); ye = nodes(elem_nodes,2);
    % Compute strain at element center (xi=0, eta=0)
    dN_dxi = 0.25 * [-1 1 1 -1];
    dN_deta=0.25*[-1 -1 1 1];
    J = [dN_dxi; dN_deta]*[xe, ye];
    invJ = inv(J);
    dN = invJ * [dN_dxi; dN_deta];
    B = zeros(3,8);
    for i_node=1:4
        B(1,2*i_node-1)=dN(1,i_node);
        B(2,2*i_node)  =dN(2,i_node);
        B(3,2*i_node-1)=dN(2,i_node);
        B(3,2*i_node)  =dN(1,i_node);
    end
    dof = reshape([2*elem_nodes-1;2*elem_nodes],1,[]);
    Ue = U(dof);
    strain(:,e) = B*Ue;
    stress(:,e) = C*strain(:,e);
end
end
