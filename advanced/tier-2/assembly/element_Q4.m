function Ke = element_Q4(coords, mat, t)

gp = [-1 1]/sqrt(3);
Ke = zeros(8,8);

for i = 1:2
    for j = 1:2
        xi  = gp(i);
        eta = gp(j);

        dNdxi = 0.25 * [ ...
            -(1-eta)  (1-eta)  (1+eta) -(1+eta)];
        dNdeta = 0.25 * [ ...
            -(1-xi)  -(1+xi)  (1+xi)   (1-xi)];

        J = [dNdxi; dNdeta] * coords;
        detJ = det(J);
        invJ = inv(J);

        dN = invJ * [dNdxi; dNdeta];

        B = zeros(3,8);
        for k = 1:4
            B(:,2*k-1:2*k) = ...
                [dN(1,k) 0;
                 0 dN(2,k);
                 dN(2,k) dN(1,k)];
        end

        Ke = Ke + B' * mat.D * B * detJ * t;
    end
end
end
