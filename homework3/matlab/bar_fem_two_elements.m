% 1D Bar Finite Element Assembly (Two Linear Elements)
% ----------------------------------------------------


clear; clc;

%% --- Problem data ---


A = 0.5;                % m^2
E = 200e9;              % Pa (200 GPa)
b = 10000.0;            % N/m (body force)
L_e = 1.0;              % element length (m)
AE = A * E;

% Two linear elements, nodes at x = 0, 1, 2
n_nodes = 3;
n_elems = 2;
connectivity = [1 2; 2 3];  % element-to-node mapping (1-based indexing in MATLAB)

%% --- Element matrices ---
% Stiffness matrix for a linear 1D element: k_e = (A*E/L) * [1 -1; -1 1]
k_elem = (AE / L_e) * [1 -1; -1 1];

% Consistent nodal load vector: f_e = b * L / 2 * [1; 1]
f_elem = b * L_e / 2 * [1; 1];

%% --- Global assembly ---
K = zeros(n_nodes, n_nodes);
F = zeros(n_nodes, 1);

for e = 1:n_elems
    conn = connectivity(e, :);
    % Assemble stiffness
    for i_local = 1:2
        i_global = conn(i_local);
        for j_local = 1:2
            j_global = conn(j_local);
            K(i_global, j_global) = K(i_global, j_global) + k_elem(i_local, j_local);
        end
    end
    % Assemble force
    for i_local = 1:2
        i_global = conn(i_local);
        F(i_global) = F(i_global) + f_elem(i_local);
    end
end

%% --- Apply Dirichlet BC: u(0) = 0 (node 1 fixed) ---
fixed_dofs = 1;
free_dofs = setdiff(1:n_nodes, fixed_dofs);

K_ff = K(free_dofs, free_dofs);
K_fi = K(free_dofs, fixed_dofs);
F_f  = F(free_dofs);

% Prescribed displacements
u = zeros(n_nodes, 1);
u(fixed_dofs) = 0.0;

% Reduced system: K_ff * u_f = F_f - K_fi * u_i
RHS = F_f - K_fi * u(fixed_dofs);

% Solve for free displacements
u_free = K_ff \ RHS;
u(free_dofs) = u_free;

%% --- Compute reactions ---
reactions = K * u - F;

%% --- Display results ---
fprintf('\nGlobal stiffness matrix K (N/m):\n');
disp(K);

fprintf('Global force vector F (N):\n');
disp(F);

fprintf('Displacements u (m) at nodes [x=0, x=1, x=2]:\n');
disp(u);

fprintf('Reactions R (N) at nodes [x=0, x=1, x=2]:\n');
disp(reactions);

fprintf('u(x=1) = %.6e m\n', u(2));
fprintf('u(x=2) = %.6e m\n', u(3));
