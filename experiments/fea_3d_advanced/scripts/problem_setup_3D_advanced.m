%% 3D FEA: Multi-material + Transient Thermal + Adaptive
Lx = 1; Ly = 0.5; Lz = 0.2;   % Geometry [m]
nx = 5; ny = 5; nz = 5;       % Initial coarse mesh

% Material properties per element
% Example: material 1 = steel, 2 = aluminum
materials = struct('E',[210e9,70e9],'nu',[0.3,0.33],'alpha',[1.2e-5,2.3e-5]);

% Load and thermal
F = 1000;                     % N
Tmax = 50;                     % Max temperature gradient
nt = 50; dt = 1;               % 50 timesteps, 1s each

% Adaptive refinement
refine_threshold = 0.8;       % % of max stress
max_iterations = 2;

% Visualization
scale_factor = 0.01;           % deformation exaggeration
results_dir = 'results/';
if ~exist(results_dir,'dir'); mkdir(results_dir); end
