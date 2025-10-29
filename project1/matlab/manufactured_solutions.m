% ------------------------------------------------------------
% Manufactured Solution Verification (1D Tapered Bar)
% ------------------------------------------------------------
clear; clc; close all;

% ------------------------------------------------------------
% Part 1: Manufactured solution and body force
% ------------------------------------------------------------
alpha = 1e-5;

u_exact = @(x) alpha * x.^2;
du_exact = @(x) 2 * alpha * x;

% --- Define bar geometry and material (consistent with tapered_bar_fem) ---
L = 2.0;    % m
E = 210e9;  % Pa
A0 = 0.02;  % m^2
A1 = 0.005; % m^2
A_func = @(x) A0 + (A1 - A0) .* (x ./ L);

% --- Package inputs as struct ---
inputs.L = L;
inputs.E = E;
inputs.A = A_func;

% --- Define body force: b(x) = -(1/A)*d/dx(E*A*du/dx) ---
body_force = @(x, inputs) ...
    -(1 ./ inputs.A(x)) .* (inputs.E .* ...
    (((inputs.A(x+1e-6) - inputs.A(x-1e-6)) ./ (2e-6)) .* du_exact(x) + ...
     inputs.A(x) .* (2*alpha)));

% ------------------------------------------------------------
% Part 2: Finite Element Solver with body force
% ------------------------------------------------------------
function [nodes, u] = fem_solve(inputs, n_el, body_force, alpha)
    n_nodes = n_el + 1;
    nodes = linspace(0, inputs.L, n_nodes);
    elems = [(1:n_el)', (2:n_el+1)'];

    K = zeros(n_nodes);
    F = zeros(n_nodes, 1);

    % Gaussian quadrature (2-point)
    xi_q = [-1/sqrt(3), 1/sqrt(3)];
    w_q  = [1, 1];

    for e = 1:n_el
        i = elems(e,1);
        j = elems(e,2);
        x1 = nodes(i);
        x2 = nodes(j);
        le = x2 - x1;

        % Element stiffness matrix
        Ae = (inputs.A(x1) + inputs.A(x2)) / 2;
        ke = inputs.E * Ae / le * [1, -1; -1, 1];
        K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + ke;

        % Element load vector (body force)
        fe = zeros(2,1);
        for q = 1:2
            xi = xi_q(q);
            w  = w_q(q);
            % Map xi to x
            x = 0.5*((1 - xi)*x1 + (1 + xi)*x2);
            N = [(1 - xi)/2, (1 + xi)/2];
            fe = fe + N' * inputs.A(x) * body_force(x, inputs) * (le/2) * w;
        end
        F(i:i+2-1) = F(i:i+2-1) + fe;
    end

    % Boundary conditions: u(0)=0
    u = zeros(n_nodes,1);
    fixed_dofs = 1;
    free_dofs = 2:n_nodes;

    % Reduced system
    K_ff = K(free_dofs, free_dofs);
    F_f  = F(free_dofs);

    % Solve
    u_f = K_ff \ F_f;
    u(free_dofs) = u_f;
end

% ------------------------------------------------------------
% Part 3: Run and compare
% ------------------------------------------------------------
n_el = 10;
[nodes, u_num] = fem_solve(inputs, n_el, body_force, alpha);

u_ex = u_exact(nodes);

% Relative L2 error
error = norm(u_num - u_ex) / norm(u_ex);
fprintf('Relative L2 error = %.6e\n', error);

% ------------------------------------------------------------
% Part 4: Plot comparison
% ------------------------------------------------------------
figure('Position',[100 100 600 400]);
plot(nodes, u_ex*1e3, 'k--', 'LineWidth',1.2); hold on;
plot(nodes, u_num*1e3, 'ro-', 'LineWidth',1.3, 'MarkerFaceColor','r');
xlabel('x (m)');
ylabel('u (mm)');
title('MMS Verification of FEM');
legend('Exact (MMS)', 'FEM', 'Location', 'NorthWest');
grid on;
