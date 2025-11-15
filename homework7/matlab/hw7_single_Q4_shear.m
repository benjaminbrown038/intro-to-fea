function hw7_single_Q4_shear
% ------------------------------------------------------------
% ENGR 4350/6350 — Homework 7
% Single Q4 plane-stress element
% Left edge fixed, uniform shear traction tau applied on right edge
% ------------------------------------------------------------

%% Inputs
E   = 200e9;      % Young's modulus [Pa]
nu  = 0.30;       % Poisson ratio
L   = 5.0;        % Length [m]
h   = 1.0;        % Height [m]
b   = 1.0;        % Thickness [m]
tau = 200e3;      % Applied shear traction σ_xy [Pa]

%% Plane stress constitutive matrix
D = (E/(1-nu^2))*[ 1,   nu,          0;
                   nu,  1,           0;
                   0,   0,  (1-nu)/2 ];

%% Node coordinates (Q4)
X = [ 0, 0;
      L, 0;
      L, h;
      0, h ];

ndof = 8;         
K = zeros(ndof);
f = zeros(ndof,1);

%% 2×2 Gauss quadrature
gp = [-1/sqrt(3), 1/sqrt(3)];
w  = [1, 1];

%% ---------------- ELEMENT STIFFNESS ------------------------------------
for i = 1:2
    xi = gp(i); wi = w(i);
    for j = 1:2
        eta = gp(j); wj = w(j);

        [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

        % Jacobian
        J = [dN_dxi dN_deta]' * X;
        detJ = det(J);
        if detJ <= 0
            error('Non-positive detJ — check element orientation.');
        end
        invJ = inv(J);

        % Gradients wrt x,y
        dN_parent = [dN_dxi, dN_deta];
        dN = dN_parent * invJ';

        % B matrix
        B = zeros(3,8);
        for a = 1:4
            B(1, 2*a-1) = dN(a,1);
            B(2, 2*a  ) = dN(a,2);
            B(3, 2*a-1) = dN(a,2);
            B(3, 2*a  ) = dN(a,1);
        end

        % Stiffness accumulation
        K = K + (B' * D * B) * b * detJ * wi * wj;
    end
end

%% ---------------- RIGHT EDGE TRACTION (xi = +1) -------------------------
tvec = [0; tau];

for j = 1:2
    eta = gp(j); wj = w(j);
    xi = +1;

    [N, ~, dN_deta] = shape_Q4(xi, eta);

    % Edge metric
    x_eta = [sum(dN_deta .* X(:,1)), sum(dN_deta .* X(:,2))];
    ds_deta = norm(x_eta);

    % Consistent nodal loads
    for a = 1:4
        f(2*a-1:2*a) = f(2*a-1:2*a) + N(a)*tvec*b*ds_deta*wj;
    end
end

%% ---------------- APPLY BOUNDARY CONDITIONS -----------------------------
fixed = [1 2 7 8];
free  = setdiff(1:ndof, fixed);

u = zeros(ndof,1);
u(free) = K(free,free) \ f(free);

%% ---------------- OUTPUT RESULTS ----------------------------------------
U = reshape(u,2,4);   % 4×2 matrix of nodal [ux, uy]

fprintf('\nNodal Displacements (meters):\n');
for a = 1:4
    fprintf('Node %d: ux = %+ .6e,   uy = %+ .6e\n', a, U(a,1), U(a,2));
end

tip_y = 0.5*(U(2,2) + U(3,2));
fprintf('\nTip (right edge) average y-deflection = %.6e m\n', tip_y);

%% ---------------- OPTIONAL DEFORMED SHAPE PLOT --------------------------
scale = 2000;
XY_def = X + scale*[U(:,1), U(:,2)];

figure; hold on; axis equal; grid on;
patch('Faces',[1 2 3 4],'Vertices',X,'FaceColor','none','LineWidth',2);
patch('Faces',[1 2 3 4],'Vertices',XY_def,'FaceColor','none','LineStyle','--','LineWidth',2);
legend('Original','Deformed (scaled)');
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Q4 Beam Under Shear (tau = %.0f kPa)', tau/1e3));

end

% ------------------------------------------------------------------------
% Q4 SHAPE FUNCTIONS
% ------------------------------------------------------------------------
function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)

N = 0.25 * [(1-xi)*(1-eta);
            (1+xi)*(1-eta);
            (1+xi)*(1+eta);
            (1-xi)*(1+eta)];

dN_dxi  = 0.25 * [ -(1-eta);
                    (1-eta);
                    (1+eta);
                   -(1+eta) ];

dN_deta = 0.25 * [ -(1-xi);
                   -(1+xi);
                    (1+xi);
                    (1-xi) ];
end
