% ========================================================================
% ENGR 4350/6350 — Homework 1
% Finite Element Shape Functions Validation (1D)
% Trey Brown — Fall 2025
% ========================================================================

% This is now a script (no function header). Just run:
% >> shape_functions_1D

xi = linspace(-1,1,200);

% --------------------------------------------------------------
% 2-node Linear Element
% --------------------------------------------------------------

% nodel values
d_lin = [0; 2];   % nodal values at xi = -1 and xi = 1

% evaluate 
for k=1:length(xi)
    N = [(1 - xi(k))/2, (1 + xi(k))/2];
    B = [-0.5, 0.5];
    u_lin(k)  = N * d_lin;
    du_lin(k) = B * d_lin;
end

figure;
subplot(1,2,1);
plot(xi, u_lin, 'b','LineWidth',1.5); hold on;
plot([-1,1], d_lin, 'ro','MarkerFaceColor','r');
title('u(x) = x + 1'); xlabel('\xi'); ylabel('u(\xi)'); grid on;
subplot(1,2,2);
plot(xi, du_lin,'r','LineWidth',1.5);
title('du/dx = 1'); xlabel('\xi'); ylabel('du/dx'); grid on;


% --------------------------------------------------------------
% 3-node Quadratic Element
% --------------------------------------------------------------

% nodal values
d_quad = [-1; -1; 1];   % nodal values at xi = -1, 0, 1

% evaluate
for k=1:length(xi)
    N = [0.5*xi(k)*(xi(k)-1), 1 - xi(k)^2, 0.5*xi(k)*(xi(k)+1)];
    B = [xi(k)-0.5, -2*xi(k), xi(k)+0.5];
    u_quad(k)  = N * d_quad;
    du_quad(k) = B * d_quad;
end


figure;
subplot(1,2,1);
plot(xi, u_quad,'b','LineWidth',1.5); hold on;
plot([-1,0,1], d_quad,'ro','MarkerFaceColor','r');
title('u(x) = x^2 + x - 1'); xlabel('\xi'); ylabel('u(\xi)'); grid on;
subplot(1,2,2);
plot(xi, du_quad,'r','LineWidth',1.5);
title('du/dx = 2x + 1'); xlabel('\xi'); ylabel('du/dx'); grid on;
