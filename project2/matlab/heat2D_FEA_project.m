function heat2D_FEA_project()
% ENGR 4350/6350 — Project 2: 2-D steady heat conduction (FEM)
% 1. Inputs 
% 2. Dirichlet BCs
% 3. Neumann BC
% 4. Mesh
% 5. Assembly 
% 6. Boundary Conditions
% 7. Post Processing
% 8. Mesh Generation
% 9. Element / Shape Function
% 10. Assembly 
% 11. Boundary Helpers
% 12. Dirichlet Solver 




% 1) Inputs
k  = 20.0;
Lx = 1.0;
Ly = 0.8;  
order = 1;            
nx = 3; 
ny = 2;       

% Dirichlet BC
T_left   = 400;       
use_Tleft = true;


% Neumann BC
q_top    = 1.0e4;     
use_qtop = true;

% Mesh
mesh = make_rect_mesh(Lx,Ly,nx,ny,order);
fprintf('Mesh: %d elements, %d nodes (order=%s)\n', mesh.ne, mesh.nn, tern(order==1,'Q4','Q9'));

% Assembly
[K,F] = assemble_global_conduction(mesh,k);

% Boundary Conditions
bc = struct;
if use_Tleft
    left_nodes = find(abs(mesh.X(:)) < 1e-12);
    bc.dir_nodes = left_nodes(:);
    bc.dir_vals  = T_left*ones(numel(left_nodes),1);
else
    bc.dir_nodes = [];
    bc.dir_vals  = [];
end


if use_qtop
    top_edge = boundary_edge(mesh,'top');
    F = apply_edge_flux(F,mesh,top_edge,@(x,y) q_top);
end

[T, free, fixed_pack] = solve_with_dirichlet(K,F,bc);

% Post-Processing
[gpXY, fluxXY, dTdx_gp, dTdy_gp] = post_flux_on_gauss(mesh,k,T);
report_and_plots(mesh,T,gpXY,fluxXY,dTdx_gp,dTdy_gp,order);
disp_table_T(mesh,T);

writematrix([(1:mesh.nn).', mesh.X(:), mesh.Y(:), T(:)], 'nodal_temperatures.csv');
fprintf('Saved nodal temperatures to nodal_temperatures.csv\n');
end

% Mesh Generation
function mesh = make_rect_mesh(Lx,Ly,nx,ny,order)
if order==1
    [Xg,Yg] = ndgrid(linspace(0,Lx,nx+1), linspace(0,Ly,ny+1));
    node_id = @(i,j) (j-1)*(nx+1) + i;
    conn = zeros(nx*ny,4); e=0;
    for j=1:ny
        for i=1:nx
            e=e+1;
            n1=node_id(i,  j);
            n2=node_id(i+1,j);
            n3=node_id(i+1,j+1);
            n4=node_id(i,  j+1);
            conn(e,:)=[n1 n2 n3 n4];
        end
    end
    X=Xg(:); 
    Y=Yg(:);
elseif order==2
    xs = linspace(0,Lx,2*nx+1);
    ys = linspace(0,Ly,2*ny+1);
    [Xg,Yg] = ndgrid(xs,ys);
    node_id = @(i,j) (j-1)*(2*nx+1)+i;
    conn = zeros(nx*ny,9); 
    e=0;
    for j=1:ny
        for i=1:nx
            e=e+1; i0=2*(i-1)+1; j0=2*(j-1)+1;
            n11=node_id(i0,  j0  ); n12=node_id(i0+1,j0  ); n13=node_id(i0+2,j0  );
            n21=node_id(i0,  j0+1); n22=node_id(i0+1,j0+1); n23=node_id(i0+2,j0+1);
            n31=node_id(i0,  j0+2); n32=node_id(i0+1,j0+2); n33=node_id(i0+2,j0+2);
            conn(e,:)=[n11 n13 n33 n31 n12 n23 n32 n21 n22];
        end
    end
    X=Xg(:); 
    Y=Yg(:);
else
    error('order must be 1 (Q4) or 2 (Q9)');
end
mesh.X=X; 
mesh.Y=Y; 
mesh.conn=conn;
mesh.nn=numel(X); 
mesh.ne=size(conn,1); 
mesh.order=order;
end

% Element / Shape Function
function [N,dN_dxi,dN_deta,wts,xis,etas]=element_shape(order,rule)
if order==1
    if nargin<2, rule=2; end
    [gp1,w1]=gauss1D(rule);
    [xis,etas]=meshgrid(gp1,gp1);
    [wx,wy]=meshgrid(w1,w1);
    wts=wx(:).*wy(:); 
    xis=xis(:); 
    etas=etas(:);
    m=numel(xis);
    N=zeros(m,4); 
    dN_dxi=zeros(m,4); 
    dN_deta=zeros(m,4);
    for g=1:m
        xi=xis(g); 
        eta=etas(g);
        N(g,:)=0.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
        dN_dxi(g,:) =0.25*[-(1-eta),(1-eta),(1+eta),-(1+eta)];
        dN_deta(g,:)=0.25*[-(1-xi),-(1+xi),(1+xi),(1-xi)];
    end
else
    if nargin<2, rule=3; end
    [gp1,w1]=gauss1D(rule);
    [xis,etas]=meshgrid(gp1,gp1);
    [wx,wy]=meshgrid(w1,w1);
    wts=wx(:).*wy(:); 
    xis=xis(:); 
    etas=etas(:);
    m=numel(xis);
    N=zeros(m,9); 
    dN_dxi=zeros(m,9); 
    dN_deta=zeros(m,9);
    for g=1:m
        xi=xis(g); 
        eta=etas(g);
        L=@(t)[0.5*t.*(t-1),(1-t.^2),0.5*t.*(t+1)];
        dL=@(t)[t-0.5,-2*t,t+0.5];
        Lx=L(xi); 
        Ly=L(eta); 
        dLx=dL(xi); 
        dLy=dL(eta);
        M=[1 5 2;8 9 6;4 7 3];
        for r=1:3
            for c=1:3
                id=M(r,c);
                N(g,id)=Lx(c)*Ly(r);
                dN_dxi(g,id)=dLx(c)*Ly(r);
                dN_deta(g,id)=Lx(c)*dLy(r);
            end
        end
    end
end
end

function [gp,w]=gauss1D(n)
switch n
    case 1, gp=0; w=2;
    case 2, gp=[-1 1]/sqrt(3); w=[1 1];
    case 3, gp=[-sqrt(3/5),0,sqrt(3/5)]; w=[5/9,8/9,5/9];
    otherwise, error('use n=1..3');
end
end

% Assembly 
function [K,F]=assemble_global_conduction(mesh,k)
nn=mesh.nn; 
ne=mesh.ne;
K=zeros(nn,nn); 
F=zeros(nn,1);
for e=1:ne
    nodes=mesh.conn(e,:); 
    xe=mesh.X(nodes); 
    ye=mesh.Y(nodes);
    if mesh.order==1, 
        rule=2; 
    else, r
        ule=3; 
    end
    [N,dN_dxi,dN_deta,wts,xis,etas]=element_shape(mesh.order,rule);
    m=numel(wts); 
    Ke=zeros(numel(nodes));
    for g=1:m
        J=[dN_dxi(g,:)*xe, dN_deta(g,:)*xe;
           dN_dxi(g,:)*ye, dN_deta(g,:)*ye];
        detJ=det(J);
        if detJ<=0, error('Non-positive detJ in element %d',e); 
        end
        invJ=inv(J);
        dNdxdy=[dN_dxi(g,:);dN_deta(g,:)]'*invJ;
        B=dNdxdy.';
        Ke=Ke+(k*(B.'*B))*detJ*wts(g);
    end
    K(nodes,nodes)=K(nodes,nodes)+Ke;
end
end

% Boundary Helpers
function edge = boundary_edge(mesh,which)
edge = [];
tol = 1e-12;
for e=1:mesh.ne
    nodes = mesh.conn(e,:);
    xe = mesh.X(nodes); 
    ye = mesh.Y(nodes);
    switch which
        case 'bottom'
            if all(abs(ye([1 2]) - 0) < tol), edge = [edge; e 1]; end
        case 'right'
            if all(abs(xe([2 3]) - max(mesh.X)) < tol), edge = [edge; e 2]; end
        case 'top'
            if all(abs(ye([3 4]) - max(mesh.Y)) < tol), edge = [edge; e 3]; end
        case 'left'
            if all(abs(xe([4 1]) - 0) < tol), edge = [edge; e 4]; end
        otherwise
            error('which ∈ {bottom,right,top,left}');
    end
end
end

function F = apply_edge_flux(F,mesh,edge_list,qn_fun)
for r=1:size(edge_list,1)
    e = edge_list(r,1); 
    ledge = edge_list(r,2);
    nodes = mesh.conn(e,:);
    xe = mesh.X(nodes); 
    ye = mesh.Y(nodes);
    [gp,w] = gauss1D(2);
    for g=1:numel(gp)
        s = gp(g); 
        wt = w(g);
        switch ledge
            case 1, N = [(1-s)/2 (1+s)/2 0 0];   
                ids=[1 2];
            case 2, N = [0 (1-s)/2 (1+s)/2 0];   
                ids=[2 3];
            case 3, N = [0 0 (1-s)/2 (1+s)/2];   
                ids=[3 4];
            case 4, N = [(1+s)/2 0 0 (1-s)/2];   
                ids=[4 1];
        end
        xg = N * xe(:); 
        yg = N * ye(:);
        qn = qn_fun(xg,yg);
        dxds = 0.5*(xe(ids(2))-xe(ids(1)));
        dyds = 0.5*(ye(ids(2))-ye(ids(1)));
        Jedge = hypot(dxds,dyds);
        F(nodes(1:4)) = F(nodes(1:4)) + N.' * (qn * Jedge * wt);
    end
end
end

% Dirichlet Solver  
function [T,free,fixed_pack]=solve_with_dirichlet(K,F,bc)
nn=size(K,1);
T=nan(nn,1);
fixed=unique(bc.dir_nodes(:));
T(fixed)=bc.dir_vals(:);
free=setdiff((1:nn)',fixed);
F_mod=F(free)-K(free,fixed)*T(fixed);
K_mod=K(free,free);
T(free)=K_mod\F_mod;
fixed_pack.fixed=fixed; 
fixed_pack.Tfixed=T(fixed);
end


function [gpXY,fluxXY,dTdx_gp,dTdy_gp]=post_flux_on_gauss(mesh,k,T)
if mesh.order==1, 
    rule=2; 
else, 
    rule=3; 
end

[~,dN_dxi,dN_deta,wts,xis,etas]=element_shape(mesh.order,rule);
gpXY=[]; 
fluxXY=[]; 
dTdx_gp=[]; 
dTdy_gp=[];
for e=1:mesh.ne
    nodes=mesh.conn(e,:); xe=mesh.X(nodes); ye=mesh.Y(nodes); Te=T(nodes);
    for g=1:numel(wts)
        J=[dN_dxi(g,:)*xe, dN_deta(g,:)*xe;dN_dxi(g,:)*ye, dN_deta(g,:)*ye];
        invJ=inv(J);
        dNdxdy=[dN_dxi(g,:);dN_deta(g,:)]'*invJ;
        gradT=dNdxdy.'*Te;
        q=-k*gradT;
        xi=xis(g); 
        eta=etas(g);
        N=0.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)].';
        xg=N.'*xe; 
        yg=N.'*ye;
        gpXY=[gpXY; xg yg]; fluxXY=[fluxXY; q.'];
        dTdx_gp=[dTdx_gp; gradT(1)]; dTdy_gp=[dTdy_gp; gradT(2)];
    end
end
end

function report_and_plots(mesh,T,gpXY,fluxXY,dTdx_gp,dTdy_gp,order)
figure('Name','Numbering'); hold on; axis equal tight;
scatter(mesh.X,mesh.Y,25,'filled'); text(mesh.X,mesh.Y,string(1:mesh.nn), 'FontSize',8);
for e=1:mesh.ne
    nodes=mesh.conn(e,:); xc=mean(mesh.X(nodes)); yc=mean(mesh.Y(nodes));
    text(xc,yc,sprintf('e%d',e),'Color',[0.85 0 0.2],'FontWeight','bold');
end


title(sprintf('Global node numbers and element labels (%s)', tern(order==1,'Q4','Q9')));
xlabel('x [m]'); ylabel('y [m]'); grid on;

figure('Name','Temperature Contour');
F=scatteredInterpolant(mesh.X,mesh.Y,T,'linear','linear');
[xg,yg]=meshgrid(linspace(min(mesh.X),max(mesh.X),101), linspace(min(mesh.Y),max(mesh.Y),101));
Tg=F(xg,yg);
contourf(xg,yg,Tg,20,'LineColor','none'); colorbar; axis equal tight;
title('Temperature field T(x,y)'); xlabel('x [m]'); ylabel('y [m]'); grid on;

figure('Name','Heat Flux');
quiver(gpXY(:,1),gpXY(:,2),fluxXY(:,1),fluxXY(:,2),'AutoScale','on');
axis equal tight; title('Heat-flux vectors q = -k ∇T (Gauss points)');
xlabel('x [m]'); ylabel('y [m]'); grid on;
end

function disp_table_T(mesh,T)
tab=table((1:mesh.nn).',mesh.X(:),mesh.Y(:),T(:), ...
    'VariableNames',{'Node','x','y','T'});
fprintf('\n=== Nodal Temperatures (K) ===\n'); disp(tab);
end

function s=tern(cond,a,b)
if cond, s=a; else, s=b; end
end
