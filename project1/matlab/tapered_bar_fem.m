% -------------------------------------------------------------------------
% ENGR 4350/6350 - Project 1
% 1D Tapered Bar Finite Element Analysis (Linear + Quadratic + MMS)
%
% Features:
%   - Linear (2-node) and Quadratic (3-node) bar elements
%   - Isoparametric shape functions N(ξ), B(ξ) = dN/dx
%   - Gauss quadrature for stiffness and body force
%   - Tapered area A(x) (linear in x)
%   - Optional Manufactured Solution (MMS) verification
%   - Tip displacement and stress evaluation
% -------------------------------------------------------------------------
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1) Problem inputs
% -------------------------------------------------------------------------
L   = 2.0;          % Bar length [m]
E   = 210e9;        % Young's modulus [Pa]
A0  = 0.02;         % Area at x = 0 [m^2]
A1  = 0.005;        % Area at x = L [m^2]

% Linear taper: A(x) = A0 + (A1 - A0)*(x/L)
A_fun = @(x) A0 + (A1 - A0)/L * x;

% Boundary conditions
u0   = 0.0;         % Displacement at x=0 (fixed)
Ptip = 1e5;         % Tip force [N] at x = L (can also use sigma * A(L))

% Element options
order   = 1;        % 1 = linear (2-node), 2 = quadratic (3-node)
n_elem  = 10;       % Number of elements (for convergence study, vary this)

% Manufactured solution toggle
use_MMS = true;     % true  = use body force from MMS, no external tip load
                    % false = physical problem: tip force Ptip, no body force

%% ------------------------------------------------------------------------
% 2) Manufactured solution setup (if used)
% -------------------------------------------------------------------------
if use_MMS
    % Choose an exact displacement field: u_exact(x) = alpha * x^2
    alpha   = 1e-5;
    u_exact = @(x) alpha * x.^2;

    % For linear A(x), derive analytic body force:
    % Strong form: -(d/dx(E A u')) = A b(x)  =>  b(x) = -(1/A)*d/dx(E A u')
    kA   = (A1 - A0)/L;
    A0f  = A0;   % just to avoid shadowing in anonymous fn
    body_force = @(x) ...
        -(E ./ (A0f + kA.*x)) .* (2*alpha*kA.*x + 2*alpha*(A0f + kA.*x));
else
    u_exact    = [];
    body_force = [];   % no body force in physical problem
end

%% ------------------------------------------------------------------------
% 3) Build mesh and solve
% -------------------------------------------------------------------------
[mesh, u] = solve_tapered_bar(L, E, A_fun, n_elem, order, u0, Ptip, ...
                              use_MMS, body_force);

x_nodes = mesh.X;
u_nodes = u;

%% ------------------------------------------------------------------------
% 4) Post-processing: element stresses at Gauss points
%    (general: works for linear & quadratic elements)
% -------------------------------------------------------------------------
[x_gp, sigma_gp] = compute_stress_gauss(mesh, u, E, A_fun);

%% ------------------------------------------------------------------------
% 5) Manufactured solution error (if MMS enabled)
% -------------------------------------------------------------------------
if use_MMS && ~isempty(u_exact)
    u_ex_nodes = u_exact(x_nodes);
    rel_L2_err = norm(u_nodes - u_ex_nodes) / norm(u_ex_nodes);
    fprintf('MMS relative L2 error = %.6e\n', rel_L2_err);
else
    rel_L2_err = NaN;
end

%% ------------------------------------------------------------------------
% 6) Plots
% -------------------------------------------------------------------------
figure('Position',[100 100 1500 400]);

% (a) Area distribution
subplot(1,3,1);
xx = linspace(0,L,200);
plot(xx, A_fun(xx), 'LineWidth', 1.5);
xlabel('x [m]'); ylabel('A(x) [m^2]');
title('Cross-sectional Area');
grid on;

% (b) Displacement
subplot(1,3,2); hold on;
plot(x_nodes, u_nodes*1e3, 'o-','LineWidth',1.5);
if use_MMS && ~isempty(u_exact)
    plot(xx, u_exact(xx)*1e3, 'k--','LineWidth',1.3);
    legend('FEM','Exact (MMS)','Location','NorthWest');
else
    legend('FEM','Location','NorthWest');
end
xlabel('x [m]'); ylabel('u [mm]');
title('Nodal Displacements');
grid on;

% (c) Stresses at Gauss points
subplot(1,3,3);
plot(x_gp, sigma_gp*1e-6, 's-','LineWidth',1.5);
xlabel('x [m]'); ylabel('\sigma [MPa]');
title('Axial Stress at Gauss Points');
grid on;

sgtitle(sprintf('Tapered Bar FEM (order=%d, n_e=%d, MMS=%d)', ...
                order, n_elem, use_MMS));


% -------------------------------------------------------------------------
% Local Functions
% -------------------------------------------------------------------------

function [mesh, u] = solve_tapered_bar(L, E, A_fun, n_elem, order, ...
                                       u0, Ptip, use_MMS, body_force)
% SOLVE_TAPERED_BAR
%   Builds mesh, assembles global K and F via Gauss quadrature, applies BCs,
%   and solves Ku = F.

    % --- Build mesh (uniform in x)
    switch order
        case 1  % linear, 2-node elements
            n_nodes = n_elem + 1;
            X       = linspace(0,L,n_nodes).';
            conn    = [(1:n_elem)', (2:n_elem+1)'];  % ne x 2
        case 2  % quadratic, 3-node elements (uniform)
            % Each element has pattern [i, i+1, i+2] with overlap of 2 nodes
            n_nodes = 2*n_elem + 1;
            X       = linspace(0,L,n_nodes).';
            conn    = zeros(n_elem,3);
            for e = 1:n_elem
                n1 = 2*e - 1;
                n2 = 2*e;
                n3 = 2*e + 1;
                conn(e,:) = [n1, n2, n3];
            end
        otherwise
            error('order must be 1 (linear) or 2 (quadratic).');
    end

    mesh.X    = X;      % nodal coordinates
    mesh.conn = conn;   % connectivity
    mesh.ne   = n_elem;
    mesh.order= order;

    ndof = n_nodes;
    K    = zeros(ndof);
    F    = zeros(ndof,1);

    % --- Gauss data
    if order == 1
        % 2-pt Gauss is sufficient for linear
        xi_q = [-1/sqrt(3),  1/sqrt(3)];
        w_q  = [1,           1];
    else
        % 3-pt Gauss for quadratic
        xi_q = [-sqrt(3/5), 0, sqrt(3/5)];
        w_q  = [5/9, 8/9, 5/9];
    end

    % --- Element loop
    for e = 1:n_elem
        elem_nodes = conn(e,:);
        Xe         = X(elem_nodes);           % column vector

        if order == 1
            ke = zeros(2,2);
            fe = zeros(2,1);
        else
            ke = zeros(3,3);
            fe = zeros(3,1);
        end

        % Gauss integration over element
        for q = 1:length(xi_q)
            xi = xi_q(q);
            w  = w_q(q);

            [N, dN_dxi] = shape_functions_1D(order, xi);

            % Isoparametric mapping: x(ξ) = sum_i N_i(ξ) X_i
            x = N * Xe;                 % 1xN * Nx1 => scalar

            % Jacobian: dx/dξ = sum_i dN_i/dξ * X_i
            J = dN_dxi * Xe;            % 1xN * Nx1 => scalar

            % B(ξ) = dN/dx = dN/dξ * (1/J)
            dN_dx = dN_dxi / J;         % row vector
            B     = dN_dx;              % 1 x nen

            % Stiffness integrand: B^T E A(x) B |J|
            ke = ke + (B' * E * A_fun(x) * B) * abs(J) * w;

            % Body force (if MMS)
            if use_MMS && ~isempty(body_force)
                b_val = body_force(x);          % [N/m^3 or N/m? depends on definition]
                fe = fe + (N' * A_fun(x) * b_val) * abs(J) * w;
            end
        end

        % Assemble into global K and F
        idx = elem_nodes;
        K(idx,idx) = K(idx,idx) + ke;
        F(idx)     = F(idx)     + fe;
    end

    % --- External tip load (physical case only)
    if ~use_MMS && Ptip ~= 0
        F(end) = F(end) + Ptip;
    end

    % --- Apply Dirichlet BC at x=0: u(1) = u0
    u = zeros(ndof,1);
    fixed_dofs = 1;
    free_dofs  = setdiff(1:ndof, fixed_dofs);

    % Modify F for known displacement (if u0 != 0)
    % K*u = F => partition => K_ff u_f = F_f - K_fc u_c
    F_f = F(free_dofs) - K(free_dofs,fixed_dofs)*u0;
    K_ff= K(free_dofs, free_dofs);

    % Solve reduced system
    u_f = K_ff \ F_f;

    % Reconstruct full displacement vector
    u(fixed_dofs) = u0;
    u(free_dofs)  = u_f;
end


function [N, dN_dxi] = shape_functions_1D(order, xi)
% SHAPE_FUNCTIONS_1D
%   Returns shape function row vector N and derivative row dN_dxi at ξ.

    if order == 1
        % Linear 2-node element, nodes at ξ = -1, +1
        N      = [(1 - xi)/2, (1 + xi)/2];       % 1x2
        dN_dxi = [-0.5,        0.5];              % 1x2
    else
        % Quadratic 3-node element, nodes at ξ = -1, 0, +1
        N1 = 0.5*xi*(xi - 1);
        N2 = 1 - xi^2;
        N3 = 0.5*xi*(xi + 1);

        N = [N1, N2, N3];

        dN1_dxi = xi - 0.5;
        dN2_dxi = -2*xi;
        dN3_dxi = xi + 0.5;

        dN_dxi = [dN1_dxi, dN2_dxi, dN3_dxi];
    end
end


function [x_gp_all, sigma_gp_all] = compute_stress_gauss(mesh, u, E, A_fun)
% COMPUTE_STRESS_GAUSS
%   Evaluates axial stress σ(x) = E * B(ξ)*u_e at Gauss points in each element.

    X    = mesh.X;
    conn = mesh.conn;
    ne   = mesh.ne;
    ord  = mesh.order;

    if ord == 1
        xi_q = [-1/sqrt(3),  1/sqrt(3)];
        w_q  = [1,           1]; %#ok<NASGU>
    else
        xi_q = [-sqrt(3/5), 0, sqrt(3/5)];
        w_q  = [5/9, 8/9, 5/9]; %#ok<NASGU>
    end

    x_gp_all     = [];
    sigma_gp_all = [];

    for e = 1:ne
        elem_nodes = conn(e,:);
        Xe         = X(elem_nodes);
        ue         = u(elem_nodes);

        for q = 1:length(xi_q)
            xi = xi_q(q);

            [N, dN_dxi] = shape_functions_1D(ord, xi);

            x = N * Xe;          % physical coordinate
            J = dN_dxi * Xe;     % dx/dξ
            dN_dx = dN_dxi / J;
            B     = dN_dx;       % 1xnen

            strain = B * ue;     % scalar
            sigma  = E * strain; % scalar

            x_gp_all     = [x_gp_all; x];      %#ok<AGROW>
            sigma_gp_all = [sigma_gp_all; sigma]; %#ok<AGROW>
        end
    end
end
