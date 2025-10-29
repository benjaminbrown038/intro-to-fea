% -------------------------------------------------------------------------
% ENGR 4350/6350 - Homework 5
% Compute element stiffness [Ke] and heat generation flux {fe}
% for one 4-node quadrilateral element using 2x2 Gauss quadrature
% -------------------------------------------------------------------------

clear; clc;

% --- Given parameters ---
k = 5;          % conductivity
s = 10;         % heat generation
xy = [0 0; 0.5 0; 0.5 0.8; 0 0.5];   % nodal coordinates [x_i, y_i]

% --- Gauss points and weights (2x2 rule) ---
gp = [-1 1]/sqrt(3);
w  = [1 1];

% Initialize results
Ke = zeros(4,4);
fe = zeros(4,1);

% --- Loop over Gauss points ---
for i = 1:2
    for j = 1:2
        xi  = gp(i);
        eta = gp(j);
        wi  = w(i);
        wj  = w(j);

        % --- Shape functions (Q4) ---
        N = 0.25 * [(1 - xi)*(1 - eta);
                    (1 + xi)*(1 - eta);
                    (1 + xi)*(1 + eta);
                    (1 - xi)*(1 + eta)];

        % --- Derivatives wrt parent coords ---
        dN_dxi = 0.25 * [ -(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta)];
        dN_deta= 0.25 * [ -(1 - xi),  -(1 + xi),   (1 + xi),   (1 - xi)];

        % --- Jacobian matrix ---
        J = [dN_dxi; dN_deta] * xy;
        detJ = det(J);
        invJ = inv(J);

        % --- Derivatives wrt global coords ---
        dN = invJ * [dN_dxi; dN_deta];
        dNdx = dN(1,:);
        dNdy = dN(2,:);

        % --- B matrix (2x4) ---
        B = [dNdx; dNdy];

        % --- Element stiffness contribution ---
        Ke = Ke + (B' * k * B) * detJ * wi * wj;

        % --- Flux vector due to heat generation ---
        fe = fe + (N * s) * detJ * wi * wj;
    end
end

% --- Display results ---
disp('Element stiffness matrix [Ke] (4x4):');
disp(Ke);

disp('Flux vector due to heat generation {fe}:');
disp(fe);


