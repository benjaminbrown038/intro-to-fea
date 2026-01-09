% Generate mesh
[nodes, elements] = generate_mesh_1D(L, n);

% Assemble stiffness
K = assemble_stiffness_1D(E, A, nodes, elements);

% Force vector
Fvec = zeros(n+1,1);
Fvec(end) = F;   % Force at last node

% Apply boundary condition (u=0 at first node)
K(1,:) = 0; K(:,1) = 0; K(1,1) = 1; 
Fvec(1) = 0;

% Solve for displacements
U = K\Fvec;

% Compute strains and stresses
strain = diff(U) ./ diff(nodes);
stress = E * strain;
