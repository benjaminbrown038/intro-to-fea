function GN = GN2D(xi, eta, nn)
% GN2D  Return derivatives of shape functions wrt (xi, eta)
% Output: 2×nn  matrix  [dN/dxi; dN/deta]

switch nn
  case 4
    dN_dxi  = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta)];
    dN_deta = 0.25 * [-(1-xi), -(1+xi), (1+xi), (1-xi)];
    GN = [dN_dxi; dN_deta];

  case 9
    % 1D bases and derivatives
    Lxi   = [0.5*xi*(xi-1), 1 - xi^2, 0.5*xi*(xi+1)];
    dLxi  = [xi - 0.5, -2*xi, xi + 0.5];
    Leta  = [0.5*eta*(eta-1), 1 - eta^2, 0.5*eta*(eta+1)];
    dLeta = [eta - 0.5, -2*eta, eta + 0.5];

    dN_dxi  = [ dLxi(1)*Leta(1), dLxi(3)*Leta(1), dLxi(3)*Leta(3), dLxi(1)*Leta(3),... 
                dLxi(2)*Leta(1), dLxi(3)*Leta(2), dLxi(2)*Leta(3), dLxi(1)*Leta(2),... 
                dLxi(2)*Leta(2) ];

    dN_deta = [ Lxi(1)*dLeta(1), Lxi(3)*dLeta(1), Lxi(3)*dLeta(3), Lxi(1)*dLeta(3),... 
                Lxi(2)*dLeta(1), Lxi(3)*dLeta(2), Lxi(2)*dLeta(3), Lxi(1)*dLeta(2),...
                Lxi(2)*dLeta(2) ];

    GN = [dN_dxi; dN_deta];

  otherwise
    error('GN2D: nn must be 4 or 9');
end
end
