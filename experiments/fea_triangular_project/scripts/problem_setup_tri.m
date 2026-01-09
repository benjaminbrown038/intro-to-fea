%% 2D Triangular FEA Problem Setup
% Geometry
Lx = 1; Ly = 0.5;          % Plate dimensions [m]

% Material (isotropic, linear elastic)
E = 210e9;                 % Youngs modulus [Pa]
nu = 0.3;                  % Poissons ratio

% Load
F = 1000;                  % Force magnitude [N]

% Mesh density
nx = 20; ny = 10;          % Number of divisions along x and y

% Thermal (optional)
alpha = 0;                 % Thermal expansion coefficient
DeltaT = 0;                % Temperature change
