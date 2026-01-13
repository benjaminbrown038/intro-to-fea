function mat = material()
mat.E  = 210e9;    % Pa
mat.nu = 0.3;

% Plane stress constitutive matrix
mat.D = (mat.E/(1-mat.nu^2)) * ...
        [1 mat.nu 0;
         mat.nu 1 0;
         0 0 (1-mat.nu)/2];
end
