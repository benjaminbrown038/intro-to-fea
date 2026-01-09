%% 2D Problem Setup
% Geometry
Lx = 1; Ly = 0.5;         % Plate dimensions [m]

% Material (isotropic linear elasticity)
E = 210e9;                % Youngs modulus [Pa]
nu = 0.3;                 % Poissons ratio

% Load
F = 1000;                 % Force magnitude [N]
load_nodes = [];           % Node numbers where force applied (fill later)

% Mesh
nx = 10; ny = 5;          % Elements along x and y

% Optional: thermal expansion
alpha = 0;                % 1/K
DeltaT = 0;               % Temperature change
