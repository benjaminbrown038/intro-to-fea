function hw7_single_Q4_shear
% ------------------------------------------------------------
% ENGR 4350/6350 — Homework 7
% Single Q4 plane-stress element
% Left edge fixed, uniform shear traction tau applied on right edge
% ------------------------------------------------------------

%% Inputs
E   = 200e9;      % Youngs modulus [Pa]
nu  = 0.30;       % Poisson ratio
L   = 5.0;        % Length [m]
h   = 1.0;        % Height [m]
b   = 1.0;        % Thickness [m]
tau = 200e3;      % Applied shear traction σ_xy [Pa]

%% Plane stress constitutive matrix
D = (E/(1-nu^2))*[ 1,   nu,          0;
                   nu,  1,           0;
                   0,   0,  (1-nu)/2 ];

%% Node coordinates (Q4 in counterclockwise order)
% 1:(0,0), 2:(L,0), 3:(L,h), 4:(0,h)
X = [ 0, 0;
      L, 0;
      L, h;
      0, h ];

ndof = 8;
K = zeros(ndof, ndof);
f = zeros(ndof, 1);

%% 2x2 Gauss points and weights
gp = [-1/sqrt(3), 1/sqrt(3)];
w  = [1, 1];

% ------------------------------------------------------------
% ELEMENT STIFFNESS
% ------------------------------------------------------------
for i = 1:2
    xi = gp(i); wi = w(i);
    for j = 1:2
        eta = gp(j); wj = w(j);

        [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

        % Parent-space gradients (4x2): [dN/dxi, dN/deta]
        dNparent = [dN_dxi, dN_deta];     % 4x2

        % Jacobian J = (dNparent)^T * X
        J = (dNparent) * X;            % 2x2
        detJ = det(J);
        if detJ <= 0
            error('Non-positive detJ');
        end
        invJ = inv(J);

        % Physical gradients dN/dx, dN/dy: (4x2)
        dN = dNparent * (invJ);         % 4x2, columns = [d/dx, d/dy]

        % Build B matrix (3x8)
        B = zeros(3, ndof);
        for a = 1:4
            B(1, 2*a-1) = dN(a,1);   % dN_a/dx on ux
            B(2, 2*a)   = dN(a,2);   % dN_a/dy on uy
            B(3, 2*a-1) = dN(a,2);   % dN_a/dy on ux
            B(3, 2*a)   = dN(a,1);   % dN_a/dx on uy
        end

        % Stiffness accumulation
        K = K + (B. * D * B) * b * detJ * wi * wj;
    end
end

% ------------------------------------------------------------
% RIGHT EDGE TRACTION (xi = +1)
% ------------------------------------------------------------
tvec = [0; tau];   % traction vector [tx; ty]

for j = 1:2
    eta = gp(j); wj = w(j);
    xi = 1.0;

    [N, ~, dN_deta] = shape_Q4(xi, eta);   % N is 4x1, dN_deta is 4x1

    % Edge metric: x_eta = sum_a dN_a/deta * x_a, y_eta similarly
    x_eta = sum(dN_deta .* X(:,1));
    y_eta = sum(dN_deta .* X(:,2));
    ds_deta = norm([x_eta; y_eta]);

    % Consistent load vector assembly
    for a = 1:4
        idx = (2*a-1):(2*a);      % [ux_a, uy_a]
        f(idx) = f(idx) + N(a) * tvec * b * ds_deta * wj;
    end
end

% ------------------------------------------------------------
% APPLY BOUNDARY CONDITIONS (left edge fixed: nodes 1 and 4)
% ------------------------------------------------------------
fixed = [1, 2, 7, 8];                       % MATLAB 1-based indexing
all_dofs = 1:ndof;
free = setdiff(all_dofs, fixed);

u = zeros(ndof,1);
Kff = K(free, free);
ff  = f(free);

u(free) = Kff \ ff;

% Reshape to nodal 4x2 matrix [ux, uy]
U = reshape(u, 2, 4).;                     % rows: nodes 1..4

% Prints
fprintf('\nNodal Displacements (meters):\n');
for k = 1:4
    fprintf('Node %d: ux = %+ .6e,   uy = %+ .6e\n', k, U(k,1), U(k,2));
end

tip_y = 0.5*(U(2,2) + U(3,2));
fprintf('\nTip deflection (avg of nodes 2 & 3): %.6e m\n', tip_y);

% ------------------------------------------------------------
% PLOT DEFORMED SHAPE
% ------------------------------------------------------------
scale = 2000;                  % visualization scale
XY_def = X + scale*U;

% close polygon ordering
ord = [1 2 3 4 1];

figure; hold on; grid on; axis equal;
plot(X(ord,1),      X(ord,2),      'k-', 'LineWidth', 1.5);
plot(XY_def(ord,1), XY_def(ord,2), 'r--','LineWidth', 1.5);
legend('Original','Deformed (scaled)');
title(sprintf('Q4 Beam Under Shear (\\tau = %.0f kPa)', tau/1000));
xlabel('x [m]'); ylabel('y [m]');

end

% ============================================================
% Q4 shape functions and derivatives at (xi, eta)
% ============================================================
function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
% Shape functions (column vector 4x1)
N = 0.25 * [ (1 - xi)*(1 - eta);
             (1 + xi)*(1 - eta);
             (1 + xi)*(1 + eta);
             (1 - xi)*(1 + eta) ];

% Derivatives wrt xi (4x1)
dN_dxi = 0.25 * [ -(1 - eta);
                   (1 - eta);
                   (1 + eta);
                  -(1 + eta) ];

% Derivatives wrt eta (4x1)
dN_deta = 0.25 * [ -(1 - xi);
                   -(1 + xi);
                    (1 + xi);
                    (1 - xi) ];
end
