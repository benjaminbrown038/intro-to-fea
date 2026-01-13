function U = solve_system(K, F, fixed_dofs)

free_dofs = setdiff(1:length(F), fixed_dofs);

U = zeros(length(F),1);
U(free_dofs) = K(free_dofs,free_dofs) \ F(free_dofs);
end
