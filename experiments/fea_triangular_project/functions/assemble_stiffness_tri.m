function K = assemble_stiffness_tri(E, nu, nodes, elements)
% Assemble global stiffness for 3-node linear triangle (plane stress)
n_nodes = size(nodes,1);
K = zeros(2*n_nodes);

% Plane stress constitutive matrix
C = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

for e = 1:size(elements,1)
    elem_nodes = elements(e,:);
    xe = nodes(elem_nodes,1); ye = nodes(elem_nodes,2);

    % Compute area
    Area = polyarea(xe, ye);

    % B matrix (linear triangle)
    x1=xe(1); y1=ye(1);
    x2=xe(2); y2=ye(2);
    x3=xe(3); y3=ye(3);

    b = [y2 - y3; y3 - y1; y1 - y2];
    c = [x3 - x2; x1 - x3; x2 - x1];

    B = zeros(3,6);
    for i=1:3
        B(1,2*i-1) = b(i);
        B(2,2*i)   = c(i);
        B(3,2*i-1) = c(i);
        B(3,2*i)   = b(i);
    end
    B = B/(2*Area);

    ke = B * C * B * Area;

    % Global DOF mapping
    dof = reshape([2*elem_nodes-1; 2*elem_nodes],1,[]);
    K(dof,dof) = K(dof,dof) + ke;
end
end
