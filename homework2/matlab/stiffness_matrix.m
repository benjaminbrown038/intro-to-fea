%% ------------------------------------------------------------
% 1D Bar Element — Variable Cross-Section Area
% Computes the element stiffness matrix:
%       Ke = ∫ Bᵀ E A(x) B dx
% using 2-point Gauss quadrature.
% ------------------------------------------------------------

%% Material Property
E = 1.0e5;               % Young’s modulus (Pa)

%% Element Geometry (physical coordinates)
x1 = 0.0;                % Left node position (m)
x2 = 3.0;                % Right node position (m)
J  = (x2 - x1) / 2.0;    % Jacobian for mapping [-1,1] → [x1,x2]

%% Shape Function Derivatives (in natural coordinates ξ)
% Linear 1D bar:
% N1 = (1 - ξ)/2, N2 = (1 + ξ)/2
% dN/dξ = [-1/2   1/2]
dN_dxi = [-0.5, 0.5];

%% Convert derivatives to physical coordinates
% dN/dx = (dN/dξ) * (dξ/dx) = (dN/dξ) / J
dN_dx = dN_dxi / J;

% Strain-displacement matrix B (1x2)
B = reshape(dN_dx, 1, 2);

% Precompute BᵀB (common term in the integrand)
BTB = B' * B;

%% Cross-Sectional Area Function A(x)
% Varies quadratically with x
A = @(x) 1 + 0.05 * x.^2;

%% Gauss Quadrature (2-point)
gauss_pts = [-1/sqrt(3), 1/sqrt(3)];     % Gauss points in ξ space
weights   = [1.0, 1.0];                 % Corresponding weights

%% Numerical Integration of ∫ A(x) dξ
integral = 0.0;

for i = 1:length(gauss_pts)
    xi = gauss_pts(i);          % Gauss point (in ξ)
    w  = weights(i);            % Weight

    % Map ξ → x:   x = x(ξ)
    x = J * xi + (x1 + x2)/2;

    % Accumulate ∑ w A(x(ξ))
    integral = integral + w * A(x);
end

% Multiply by Jacobian because dx = J dξ
integral = integral * J;

%% Element Stiffness Matrix
% Ke = E * ∫ (Bᵀ B A(x)) dx
Ke = E * BTB * integral;

%% Output
disp('Element stiffness matrix Ke:');
disp(Ke);
