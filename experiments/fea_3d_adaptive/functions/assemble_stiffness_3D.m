function K = assemble_stiffness_3D(E, nu, nodes, elements)
n_nodes = size(nodes,1);
K = zeros(3*n_nodes);

C = E/((1+nu)*(1-2*nu)) * ...
    [1-nu nu nu 0 0 0;
     nu 1-nu nu 0 0 0;
     nu nu 1-nu 0 0 0;
     0 0 0 (1-2*nu)/2 0 0;
     0 0 0 0 (1-2*nu)/2 0;
     0 0 0 0 0 (1-2*nu)/2];

for e = 1:size(elements,1)
    elem_nodes = elements(e,:);
    xe = nodes(elem_nodes,1); ye = nodes(elem_nodes,2); ze = nodes(elem_nodes,3);

    V = abs(det([1 1 1 1; xe; ye; ze]))/6;

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

    ke = B * C * B * V;
    dof = reshape([3*elem_nodes-2;3*elem_nodes-1;3*elem_nodes],1,[]);
    K(dof,dof) = K(dof,dof) + ke;
end
end
