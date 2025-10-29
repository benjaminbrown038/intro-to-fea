function plane_stress_FEA_project()
% ENGR 4350/6350 Project 3 — 2D Linear Elasticity 

% 1. Material 
% Matrix
% 2. Geometry 
% 3. Mesh
% 3. Boundary Conditions
% 4. Apply Uniform Traction on Right Edge 
% 5. Assembly 
% 6. Dirichlet BCs
% 7. Post Process 
% 8. Displacement Table and CSV
% 9. Mesh and Shape Functions 
% 10. Assembly 
% 11. Post Process (Stress, Mises)
% 12. Plotting 
% 13. Node Numbering 
% 14. Displacement Vectors
% 15. Choose a scale that roughly fills 10% of domain width 
% 16. Von Mises Contour from Gauss Points
% 17. Boundary Helpers
 

% material 
E  = 210e9;
nu = 0.30;
planeStress = true;

% Matrix (Strain or Stress)
if planeStress
    D = E/(1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0,(1-nu)/2];
else 
    D = E/((1+nu)*(1-2*nu)) * [1-nu, nu, 0; nu, 1-nu, 0; 0,0,(1-2*nu)/2];
end

% Geometry / Mesh 
Lx = 1.0;   
Ly = 0.8;
nx = 3;     
ny = 2;     
order = 1;              

mesh = make_rect_mesh(Lx,Ly,nx,ny,order);
fprintf('Mesh: %d elements, %d nodes\n', mesh.ne, mesh.nn);

% Boundary Conditions 
left_nodes = find(abs(mesh.X) < 1e-12);
bc.fixed_dofs = [2*left_nodes-1; 2*left_nodes];  
bc.fixed_dofs = bc.fixed_dofs(:);
bc.fixed_vals = zeros(numel(bc.fixed_dofs),1);

% Apply Uniform Traction on Right Edge 
right_edges = boundary_edge(mesh,'right');
traction = @(x,y) [1e6; 0];

% Assembly 
[K,F] = assemble_elastic(mesh,D,traction,right_edges);

% Apply Dirichlet BCs 
ndof = 2*mesh.nn;
free = setdiff((1:ndof)', bc.fixed_dofs);
u = zeros(ndof,1);
u(bc.fixed_dofs) = bc.fixed_vals;
F_mod = F(free) - K(free,bc.fixed_dofs) * bc.fixed_vals;
K_mod = K(free,free);
u(free) = K_mod \ F_mod;

% Post Process 
[stress, vonMises, coord_gp] = stress_post(mesh,D,u);
plot_results(mesh,u,coord_gp,stress,vonMises);

% Displacement Table and CSV
disp('Displacements at each node [ux (m), uy (m)]:');
Utab = [mesh.X(:), mesh.Y(:), u(1:2:end), u(2:2:end)];
disp(array2table(Utab,'VariableNames',{'x','y','ux','uy'}));
writematrix(Utab,'nodal_displacements.csv');
fprintf('Saved nodal displacements to nodal_displacements.csv\n');
end

% Mesh & Shape Functions
function mesh = make_rect_mesh(Lx,Ly,nx,ny,order)
if order ~= 1
    error('Only Q4 is implemented in this file (set order = 1).');
end

[Xg,Yg] = ndgrid(linspace(0,Lx,nx+1), linspace(0,Ly,ny+1));
X = Xg(:);  
Y = Yg(:);
conn = zeros(nx*ny, 4); 
e = 0;

for j = 1:ny
    for i = 1:nx
        e  = e + 1;
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        conn(e,:) = [n1 n2 n3 n4];  
    end
end

mesh.X = X; 
mesh.Y = Y; 
mesh.conn = conn;
mesh.ne = size(conn,1); 
mesh.nn = numel(X); 
mesh.order = order;
end

function [gp,w] = gauss1D(n)
switch n
    case 2
        gp = [-1 1] / sqrt(3);  
        w = [1 1];
    case 3
        gp = [-sqrt(3/5) 0 sqrt(3/5)];  
        w = [5/9 8/9 5/9];
    otherwise
        error('gauss1D: use n=2 or 3');
end
end

function [N,dN_dxi,dN_deta,wts,xis,etas] = shape_Q4()
[gp,w] = gauss1D(2);
[xis,etas] = meshgrid(gp,gp);
[wx,wy] = meshgrid(w,w);
wts = (wx(:) .* wy(:));
xis = xis(:);  
etas = etas(:);

m = numel(xis);
N = zeros(m,4); 
dN_dxi = N; 
dN_deta = N;

for g = 1:m
    xi = xis(g);  
    eta = etas(g);
    N(g,:) = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN_dxi(g,:) = 0.25 * [-(1-eta), +(1-eta), +(1+eta),-(1+eta)];
    dN_deta(g,:) = 0.25 * [-(1-xi), -(1+xi), +(1+xi),+(1-xi)];
end
end

% Assembly 
function [K,F] = assemble_elastic(mesh,D,traction,right_edges)
ndof = 2*mesh.nn;
K = zeros(ndof, ndof);
F = zeros(ndof, 1);
[N, dN_dxi, dN_deta, wts, ~, ~] = shape_Q4();

for e = 1:mesh.ne
    nodes = mesh.conn(e,:);
    xe = mesh.X(nodes);  
    ye = mesh.Y(nodes);
    Ke = zeros(8,8);

    for g = 1:numel(wts)
        J = [ dN_dxi(g,:)*xe, dN_deta(g,:)*xe; dN_dxi(g,:)*ye, dN_deta(g,:)*ye];
        detJ = det(J);
        if detJ <= 0
            error('Non-positive detJ in element %d (check node ordering/mesh)', e);
        end

        invJ  = inv(J);
        dNdxdy = [dN_dxi(g,:); dN_deta(g,:)]' * invJ;
        B = zeros(3,8);

        for i = 1:4
            B(:, [2*i-1, 2*i]) = [dNdxdy(i,1)0;0dNdxdy(i,2);dNdxdy(i,2)  dNdxdy(i,1) ];
        end

        Ke = Ke + (B' * D * B) * (detJ * wts(g));

    end
    dofs = reshape([2*nodes-1; 2*nodes], 1, []);
    K(dofs, dofs) = K(dofs, dofs) + Ke;
end


for r = 1:size(right_edges,1)
    e      = right_edges(r,1);
    nodes2 = mesh.conn(e, [2 3]);   
    xe     = mesh.X(nodes2);  
    ye = mesh.Y(nodes2);
    [gp1, w1] = gauss1D(2);
    Ledge = hypot(diff(xe), diff(ye));  
    for g = 1:numel(gp1)
        s  = gp1(g);  
        wt = w1(g);
        N1 = (1 - s)/2;  
        N2 = (1 + s)/2;
        xg = N1*xe(1) + N2*xe(2);
        yg = N1*ye(1) + N2*ye(2);
        t  = traction(xg, yg);         
        Fe_local = [N1*t;N2*t ] * (Ledge/2 * wt);  
        dofs = [2*nodes2(1)-1, 2*nodes2(1), 2*nodes2(2)-1, 2*nodes2(2) ];
        F(dofs) = F(dofs) + Fe_local(:);
    end
end
end

% Post Process (Stress, Mises) 
function [stress, vonMises, coord_gp] = stress_post(mesh,D,u)
[N, dN_dxi, dN_deta, wts, xis, etas] = shape_Q4();
stress   = [];
vonMises = [];
coord_gp = [];

for e = 1:mesh.ne
    nodes = mesh.conn(e,:);
    xe = mesh.X(nodes);  ye = mesh.Y(nodes);
    ue = u(reshape([2*nodes-1; 2*nodes], 1, []).');
    for g = 1:numel(wts)
        J = [dN_dxi(g,:)*xe, dN_deta(g,:)*xe; dN_dxi(g,:)*ye, dN_deta(g,:)*ye];
        invJ   = inv(J);
        dNdxdy = [dN_dxi(g,:); dN_deta(g,:)]' * invJ;
        B = zeros(3,8);
        for i = 1:4
            B(:, [2*i-1, 2*i]) = [dNdxdy(i,1)0;0dNdxdy(i,2);dNdxdy(i,2)dNdxdy(i,1) ];
        end
        sig = D * (B * ue);   
        vm  = sqrt(sig(1)^2 - sig(1)*sig(2) + sig(2)^2 + 3*sig(3)^2);
        xi = xis(g);  eta = etas(g);
        Ng = 0.25 * [ (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta) ]';
        xg = Ng' * xe;  yg = Ng' * ye;
        coord_gp = [coord_gp; xg yg];
        stress   = [stress;   sig'];
        vonMises = [vonMises; vm];
    end
end
end

% Plotting 
function plot_results(mesh,u,coord_gp,stress,vonMises)

% Node numbering
figure('Name','Node numbering'); hold on; axis equal; box on;
scatter(mesh.X, mesh.Y, 25, 'filled');
text(mesh.X, mesh.Y, string(1:mesh.nn), 'FontSize', 8, 'VerticalAlignment','bottom');
title('Global node numbers'); xlabel('x'); ylabel('y'); grid on;

% Displacement Vectors 
ux = u(1:2:end);  uy = u(2:2:end);

% Choose a scale that roughly fills 10% of domain width
umax = max(sqrt(ux.^2 + uy.^2) + eps);
scale = 0.1 * (max(mesh.X)-min(mesh.X)) / umax;
figure('Name','Displacement vectors');
quiver(mesh.X, mesh.Y, ux, uy, scale, 'AutoScale','off');
axis equal; box on; grid on; xlabel('x'); ylabel('y');
title('Displacement vectors (scaled)');

% Von Mises Contour from Gauss Points
F = scatteredInterpolant(coord_gp(:,1), coord_gp(:,2), vonMises, 'linear', 'linear');
[xg,yg] = meshgrid(linspace(min(mesh.X),max(mesh.X),100), linspace(min(mesh.Y),max(mesh.Y),100));
vm_grid = F(xg,yg);
figure('Name','Von Mises stress');
contourf(xg, yg, vm_grid, 20, 'LineColor','none'); colorbar;
axis equal tight; box on; xlabel('x [m]'); ylabel('y [m]');
title('Von Mises stress (Pa)');
end


% Boundary Helpers
function edge = boundary_edge(mesh,which)

edge = [];
tol  = 1e-12;
for e = 1:mesh.ne
    nodes = mesh.conn(e,:);
    xe = mesh.X(nodes);  ye = mesh.Y(nodes);
    switch lower(which)
        case 'bottom'
            if all(abs(ye([1 2]) - min(mesh.Y)) < tol), edge = [edge; e 1]; end
        case 'right'
            if all(abs(xe([2 3]) - max(mesh.X)) < tol), edge = [edge; e 2]; end
        case 'top'
            if all(abs(ye([3 4]) - max(mesh.Y)) < tol), edge = [edge; e 3]; end
        case 'left'
            if all(abs(xe([4 1]) - min(mesh.X)) < tol), edge = [edge; e 4]; end
        otherwise
            error('boundary_edge: which must be one of {bottom,right,top,left}');
    end
end
end
