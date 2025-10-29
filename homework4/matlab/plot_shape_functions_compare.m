

function plot_shape_functions_compare()
% Compare 4- and 9-node 2D shape functions visually

xiplot  = linspace(-1, 1, 25);
etaplot = linspace(-1, 1, 25);
[Xi, Eta] = meshgrid(xiplot, etaplot);

for nn = [4, 9]
    theta     = zeros(size(Xi));
    dth_dxi   = zeros(size(Xi));
    dth_deta  = zeros(size(Xi));
    d         = (1:nn)';   % arbitrary nodal values (just for plotting)

    % Loop through every (xi, eta) point
    for i = 1:numel(Xi)
        xi  = Xi(i);
        eta = Eta(i);

        N   = N2D(xi, eta, nn);
        N   = reshape(N, 1, []);   % ensure row vector
        GN  = GN2D(xi, eta, nn);

        theta(i)    = N * d;
        dth_dxi(i)  = GN(1,:) * d;
        dth_deta(i) = GN(2,:) * d;
    end

    % --- Surface plot ---
    figure
    surf(Xi, Eta, theta)
    title(sprintf('Shape Function Surface (%d-node)', nn))
    xlabel('\xi'); ylabel('\eta'); zlabel('\theta')
    shading interp; grid on

    % --- Gradient field plot ---
    figure
    quiver(Xi, Eta, dth_dxi, dth_deta)
    title(sprintf('Gradient Field (%d-node)', nn))
    xlabel('\xi'); ylabel('\eta')
    axis equal; grid on
end
end
