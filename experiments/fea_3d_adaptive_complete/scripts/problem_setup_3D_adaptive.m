%% 3D Adaptive Thermal-Mechanical FEA Setup
Lx = 1; Ly = 0.5; Lz = 0.2;     % Dimensions [m]

E = 210e9; nu = 0.3;            % Material properties
alpha = 1.2e-5;                 % Thermal expansion [1/K]

F = 1000;                       % Mechanical load [N]

nx = 5; ny = 5; nz = 5;         % Initial coarse mesh

Tmax = 50;                       % Max temperature for gradient

refine_threshold = 0.8;          % refine elements with stress > 80% max
max_iterations = 2;              % number of adaptive refinement loops

scale_factor = 0.01;             % deformation exaggeration for visualization
