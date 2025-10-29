% ------------------------------------------------------------
% Tapered Bar Finite Element Analysis (1D)
% ------------------------------------------------------------
clear; clc; close all;

% ----------------------------
% Part 1: Inputs
% ----------------------------

% Bar length (m)
L = 2.0; 

% Young's modulus (Pa)
E = 210e9;          

% Area at x=0 (m^2)
A0 = 0.02; 

% Area at x=L (m^2)
A1 = 0.005;   

% Fixed end
u0_prescribed = 0.0; 

% Not used
sigma_L_prescribed = []; 

% area function
linear_taper = @(A0, A1, L) @(x) A0 + (A1 - A0) .* (x ./ L);

% area function 
exponential_taper = @(A0, A1, L) @(x) A0 .* (A1./A0) .^ (x ./ L);

% Select function (linear or exponential)
A_func = linear_taper(A0, A1, L);

% --- Summary output ---
fprintf('Tapered bar inputs summary\n');
fprintf('--------------------------\n');
fprintf('L = %.3f m\n', L);
fprintf('E = %.3e Pa\n', E);
fprintf('u(0) = %.3f\n', u0_prescribed);
if isempty(sigma_L_prescribed)
    fprintf('sigma(L) = None\n');
else
    fprintf('sigma(L) = %.3e\n', sigma_L_prescribed);
end
xs = linspace(0, L, 5);
areas = A_func(xs);
fprintf('Sample A(x):\n');
for k = 1:length(xs)
    fprintf('  A(%.4g) = %.6g\n', xs(k), areas(k));
end

% ----------------------------
% Part 2: Finite Element Solver
% ----------------------------

% Number of elements
n_el = 10; 
% Number of nodes
n_nodes = n_el + 1;
% Node coordinates
nodes = linspace(0, L, n_nodes);  
% Connectivity
elems = [(1:n_el)', (2:n_el+1)']; 

% stiffness 
K = zeros(n_nodes);  

% force 
F = zeros(n_nodes, 1); 

% stiffness matrix 
for e = 1:n_el
    i = elems(e,1);
    j = elems(e,2);
    x1 = nodes(i);
    x2 = nodes(j);
    le = x2 - x1;
    Ae = (A_func(x1) + A_func(x2)) / 2;
    ke = E * Ae / le * [1, -1; -1, 1];
    K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + ke;
end

% Boundary conditions 

% Node 1 fixed
u = zeros(n_nodes, 1);

% Node 1 fixed
fixed_dofs = 1;

% free nodes 
free_dofs = 2:n_nodes;

% force 
F(end) = 1e5;                       

% --- Solve for unknown displacements ---
K_ff = K(free_dofs, free_dofs);
F_f = F(free_dofs);
u_f = K_ff \ F_f;
u(free_dofs) = u_f;

% ----------------------------
% Part 3: Stresses
% ----------------------------
stresses = zeros(n_el,1);
for e = 1:n_el
    i = elems(e,1);
    j = elems(e,2);
    le = nodes(j) - nodes(i);
    du = u(j) - u(i);
    strain = du / le;
    stresses(e) = E * strain;
end

% ----------------------------
% Part 4: Plots
% ----------------------------
figure('Position',[100 100 1400 400]);

% (a) Area distribution
subplot(1,3,1);
xs = linspace(0, L, 200);
plot(xs, A_func(xs), 'LineWidth',1.5);
title('Cross-sectional Area A(x)');
xlabel('x (m)');
ylabel('A(x) (m^2)');
grid on;

% (b) Displacements
subplot(1,3,2);
plot(nodes, u*1e3, 'o-', 'LineWidth',1.5);
title('Nodal Displacements');
xlabel('x (m)');
ylabel('u (mm)');
grid on;

% (c) Stresses
subplot(1,3,3);
x_centers = nodes(1:end-1) + diff(nodes)/2;
plot(x_centers, stresses*1e-6, 's-r', 'LineWidth',1.5);
title('Element Stresses');
xlabel('x (m)');
ylabel('Stress (MPa)');
grid on;

sgtitle('Tapered Bar Finite Element Results');
