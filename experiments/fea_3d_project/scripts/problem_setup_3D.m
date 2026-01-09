%% 3D Tetrahedral Thermal-Mechanical FEA Problem Setup
Lx = 1; Ly = 0.5; Lz = 0.2;     % Cube/plate dimensions [m]

E = 210e9; nu = 0.3;            % Material properties
alpha = 1.2e-5;                 % Thermal expansion coefficient [1/K]

F = 1000;                       % Applied mechanical load

% Mesh density (divisions along x, y, z)
nx = 5; ny = 5; nz = 5;

Tmax = 50;                       % Maximum temperature for gradient
