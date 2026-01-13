function U = solve_system(K, F, fixed_dofs)

free = setdiff(1:length(F), fixed_dofs);

U = zeros(length(F),1);
U(free) = K(free,free) \ F(free);
end
