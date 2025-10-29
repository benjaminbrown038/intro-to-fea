function N = N2D(xi, eta, nn)


switch nn
case 4
  N = 0.25 * [ (1-xi)*(1-eta), 
                 (1+xi)*(1-eta), 
                 (1+xi)*(1+eta), 
                 (1-xi)*(1+eta) ];

case 9
  Lxi  = [0.5*xi*(xi-1), 1 - xi^2, 0.5*xi*(xi+1)];
  Leta = [0.5*eta*(eta-1), 1 - eta^2, 0.5*eta*(eta+1)];

  N = [ Lxi(1)*Leta(1), Lxi(3)*Leta(1), Lxi(3)*Leta(3), Lxi(1)*Leta(3), ...
      Lxi(2)*Leta(1), Lxi(3)*Leta(2), Lxi(2)*Leta(3), Lxi(1)*Leta(2), ...
      Lxi(2)*Leta(2) ];
otherwise
  error('N2D: nn must be 4 or 9');
end
end
