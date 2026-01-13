function Ke = element_stiffness(E, nu, t, coords)

D = (E/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

xi = [-1 1 1 -1]/sqrt(3);
eta = [-1 -1 1 1]/sqrt(3);

Ke = zeros(8,8);

for i = 1:4
    dNdxi = 0.25*[-(1-eta(i)) (1-eta(i)) (1+eta(i)) -(1+eta(i))];
    dNdeta= 0.25*[-(1-xi(i)) -(1+xi(i)) (1+xi(i)) (1-xi(i))];

    J = [dNdxi; dNdeta]*coords;
    detJ = det(J);
    invJ = inv(J);

    dN = invJ*[dNdxi; dNdeta];

    B = zeros(3,8);
    for k = 1:4
        B(:,2*k-1:2*k) = [dN(1,k) 0; 0 dN(2,k); dN(2,k) dN(1,k)];
    end

    Ke = Ke + B*D*B*detJ*t;
end
end
