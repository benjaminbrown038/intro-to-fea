function mat = material()
mat.E  = 210e9;
mat.nu = 0.3;

mat.D = (mat.E/(1-mat.nu^2)) * ...
        [1 mat.nu 0;
         mat.nu 1 0;
         0 0 (1-mat.nu)/2];
end
