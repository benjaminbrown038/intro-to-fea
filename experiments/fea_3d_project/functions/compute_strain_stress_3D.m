function [strain, stress] = compute_strain_stress_3D(U, nodes, elements, E, nu, alpha, Tmax)
n_elem = size(elements,1);
strain = zeros(6,n_elem);
stress = zeros(6,n_elem);

C = E/((1+nu)*(1-2*nu)) * ...
    [1-nu nu nu 0 0 0;
     nu 1-nu nu 0 0 0;
     nu nu 1-nu 0 0 0;
     0 0 0 (1-2*nu)/2 0 0;
     0 0 0 0 (1-2*nu)/2 0;
     0 0 0 0 0 (1-2*nu)/2];

% Non-uniform temperature: gradient along z
z_centers = mean(nodes(elements,3),2);
DeltaT_elem = Tmax * z_centers / max(nodes(:,3));

for e = 1:n_elem
    elem_nodes = elements(e,:);
    xe = nodes(elem_nodes,1); ye = nodes(elem_nodes,2); ze = nodes(elem_nodes,3);

    % Compute tetrahedral volume and B as above
    D = [1 xe(1) ye(1) ze(1);
         1 xe(2) ye(2) ze(2);
         1 xe(3) ye(3) ze(3);
         1 xe(4) ye(4) ze(4)];
    invD = inv(D);
    a = invD(2,:); b = invD(3,:); c = invD(4,:);
    B = zeros(6,12);
    for i = 1:4
        B(:,3*i-2:3*i) = [a(i) 0 0; 0 b(i) 0; 0 0 c(i); b(i) a(i) 0; 0 c(i) b(i); c(i) 0 a(i)];
    end

    dof = reshape([3*elem_nodes-2;3*elem_nodes-1;3*elem_nodes],1,[]);
    Ue = U(dof);

    strain(:,e) = B * Ue;
    strain_thermal = alpha * DeltaT_elem(e) * [1;1;1;0;0;0];  % thermal strain
    stress(:,e) = C * (strain(:,e) - strain_thermal);
end
end
