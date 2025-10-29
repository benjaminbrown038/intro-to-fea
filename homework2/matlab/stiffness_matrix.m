
% Material and geometry parameters

% Youngs modulus
E = 1.0e5;               

% Element domain
x1 = 0.0;                

% Element domain
x2 = 3.0;                

% Jacobian
J = (x2 - x1) / 2.0;     

% Shape function derivatives 
dN_dxi = [-0.5, 0.5];

% Convert 
dN_dx = dN_dxi / J;
B = reshape(dN_dx, 1, 2);   % Row vector

% Compute 
BTB = B' * B;

% Area function
A = @(x) 1 + 0.05 * x.^2;


gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
weights = [1.0, 1.0];

% Integral 
integral = 0.0;
for i = 1:length(gauss_pts)
    xi = gauss_pts(i);
    w = weights(i);
    x = (x2 - x1)/2 * xi + (x1 + x2)/2;
    integral = integral + w * A(x);
end

integral = integral * J;  

% Stiffness matrix
Ke = E * BTB * integral;

% Display result
disp('Element stiffness matrix Ke:');
disp(Ke);
