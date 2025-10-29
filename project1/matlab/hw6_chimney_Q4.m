function hw6_chimney_Q4()
% ENGR 4350/6350 — Homework 6: Chimney Conduction (Q4 FEM)


k_brick    = 0.9;   

k_concrete = 2.0;   

T_in  = 140.0;   % inside (hot gas)

T_out = 10.0;    % outside (ambient)

% Node coordinates (approx. 1/8 circular region)
X = [0.30 0.35 0.40;
     0.30 0.40 0.50;
     0.30 0.45 0.60];

Y = [0.40 0.40 0.40;
     0.50 0.50 0.50;
     0.60 0.60 0.60];

nodes = [X(:), Y(:)];

ndof = size(nodes,1);

% Element connectivity (counterclockwise)
conn = [
    1 4 5 2;
    2 5 6 3;
    4 7 8 5;
    5 8 9 6];

% Material assignment per element
k_elem = [k_brick; k_concrete; k_brick; k_concrete];

% Boundary conditions — only inner/outer arcs fixed T
tol = 1e-12;
x = nodes(:,1); 
y = nodes(:,2);

% diagonal nodes approximate circular arcs
on_diag = abs(x - y) < tol;

inner_arc_nodes = find(on_diag & (x <= 0.40 + tol));   % inner arc (T_in)
outer_arc_nodes = find(on_diag & (x >= 0.60 - tol));   % outer arc (T_out)

BC_nodes = [inner_arc_nodes; outer_arc_nodes];
BC_vals  = [repmat(T_in,  numel(inner_arc_nodes),1);
            repmat(T_out, numel(outer_arc_nodes),1)];

% Assembly
K = zeros(ndof, ndof);
[gp, gw] = gauss2x2();

for e = 1:size(conn,1)
    en = conn(e,:);          
    xe = nodes(en,1);        
    ye = nodes(en,2);        
    ke = zeros(4,4);
    kappa = k_elem(e);

    for ig = 1:4
        xi  = gp(ig,1);  
        eta = gp(ig,2);  
        w   = gw(ig);

        [N, dN_dxi, dN_deta] = N_Q4(xi, eta);

        J = [dN_dxi'*xe,  dN_dxi'*ye;
             dN_deta'*xe, dN_deta'*ye];
        detJ = det(J);
        if detJ <= 0
            error('Non-positive detJ in element %d', e);
        end
        invJ = inv(J);

        dN = [dN_dxi dN_deta] * invJ.';   % (ξ,η) → (x,y)
        B  = dN.';                         
        ke = ke + (B.' * kappa * eye(2) * B) * detJ * w;
    end

    K(en,en) = K(en,en) + ke;
end

% Apply BCs and solve
f = zeros(ndof,1);
free = true(ndof,1);
free(BC_nodes) = false;

Kff = K(free,free);
Kfc = K(free,~free);

% reorder BC_vals to correspond to ~free
uc = nan(sum(~free),1);
[~,loc] = ismember(BC_nodes, find(~free));
uc(loc) = BC_vals;

rhs = f(free) - Kfc * uc;
uf = Kff \ rhs;

u = zeros(ndof,1);
u(free) = uf;  
u(~free) = uc;

% Output
fprintf('\nNodal Temperatures (°C):\n');
for n = 1:ndof
    fprintf('  Node %d (%.3f, %.3f):  T = %8.4f\n', n, nodes(n,1), nodes(n,2), u(n));
end

fprintf('\nNode numbering (x, y):\n');
for n = 1:ndof
    fprintf('  %2d: (%.3f, %.3f)\n', n, nodes(n,1), nodes(n,2));
end
fprintf('\nElement connectivity (Q4):\n');
for e = 1:size(conn,1)
    fprintf('  e%-2d: [%d %d %d %d]\n', e, conn(e,1),conn(e,2),conn(e,3),conn(e,4));
end

% 7) Visualization
figure('Name','HW6 Chimney — Temperature Field');
trisurf_draped_as_contour(nodes, conn, u);
axis equal tight; view(2); colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Temperature Distribution (°C) — Q4 FEM');

% Node and element numbering
figure('Name','Node & Element Numbering'); hold on; axis equal tight;
patch('Faces',conn,'Vertices',nodes,'FaceColor','none','EdgeColor',[0.4 0.4 0.4]);
for n = 1:ndof
    text(nodes(n,1), nodes(n,2), sprintf('%d',n), 'Color','b', 'FontSize',10, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
for e = 1:size(conn,1)
    c = mean(nodes(conn(e,:),:),1);
    text(c(1), c(2), sprintf('e%d',e), 'Color','r', 'FontSize',10, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
xlabel('x [m]'); ylabel('y [m]');
title('Node (blue) and Element (red) Numbering');
grid on;

end

% Helper Functions
function [N, dN_dxi, dN_deta] = N_Q4(xi, eta)
N1 = 0.25*(1 - xi)*(1 - eta);
N2 = 0.25*(1 + xi)*(1 - eta);
N3 = 0.25*(1 + xi)*(1 + eta);
N4 = 0.25*(1 - xi)*(1 + eta);
N  = [N1 N2 N3 N4];
dN_dxi  = 0.25*[-(1 - eta); +(1 - eta); +(1 + eta); -(1 + eta)];
dN_deta = 0.25*[-(1 - xi ); -(1 + xi ); +(1 + xi ); +(1 - xi )];
end

function [gp, gw] = gauss2x2()
g = 1/sqrt(3);
gp = [-g -g; +g -g; +g +g; -g +g];
gw = [1;1;1;1];
end

function trisurf_draped_as_contour(nodes, conn, u)
tri = [conn(:,[1 2 3]); conn(:,[1 3 4])];
trisurf(tri, nodes(:,1), nodes(:,2), u, 'EdgeAlpha', 0.4);
shading interp;
end
