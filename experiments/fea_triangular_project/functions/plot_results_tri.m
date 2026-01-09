function plot_results_tri(nodes, elements, stress)
% Compute von Mises stress
sigma_x = stress(1,:);
sigma_y = stress(2,:);
tau_xy  = stress(3,:);
sigma_vm = sqrt(sigma_x.^2 - sigma_x.*sigma_y + sigma_y.^2 + 3*tau_xy.^2);

figure;
trisurf(elements, nodes(:,1), nodes(:,2), sigma_vm, 'FaceColor','interp','EdgeColor','k');
colorbar; 
title('Von Mises Stress');
xlabel('x [m]'); ylabel('y [m]');
view(2); axis equal; grid on;
end
