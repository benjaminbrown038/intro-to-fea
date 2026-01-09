function K = assemble_stiffness_2D(E, nu, nodes, elements)
% 2D linear plane stress quad elements

n_nodes = size(nodes,1);
K = zeros(2*n_nodes);

% Plane stress constitutive matrix
C = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% Loop over elements
for e = 1:size(elements,1)
    elem_nodes = elements(e,:);
    xe = nodes(elem_nodes,1); ye = nodes(elem_nodes,2);

    % 4-node bilinear quad stiffness (simplified using 2x2 Gauss)
    ke = zeros(8,8);
    gauss = [-1/sqrt(3), 1/sqrt(3)];
    for xi = gauss
        for eta = gauss
            % Shape functions derivatives
            dN_dxi = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta=0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];
            J = [dN_dxi; dN_deta]*[xe, ye];
            detJ = det(J);
            invJ = inv(J);
            dN = invJ * [dN_dxi; dN_deta];
            B = zeros(3,8);
            for i_node=1:4
                B(1,2*i_node-1)=dN(1,i_node);
                B(2,2*i_node)  =dN(2,i_node);
                B(3,2*i_node-1)=dN(2,i_node);
                B(3,2*i_node)  =dN(1,i_node);
            end
            ke = ke + B' * C * B * detJ;
        end
    end
    
    % Assemble into global K
    dof = reshape([2*elem_nodes-1;2*elem_nodes],1,[]);
    K(dof,dof) = K(dof,dof) + ke;
end
end
