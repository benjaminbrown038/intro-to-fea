function plot_results(nodes, elements, U, stress)

scale = 50;
def = nodes + scale*[U(1:2:end), U(2:2:end)];

figure;
trisurf(elements, def(:,1), def(:,2), zeros(size(def,1),1), ...
        stress(:,1), 'EdgeColor','k');
view(2); axis equal; colorbar;
title('\sigma_x (Triangular Elements)');
end
