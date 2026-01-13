function Ke = element_T3(coords, mat, t)

x = coords(:,1); y = coords(:,2);

A = polyarea(x,y);

b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];

B = 1/(2*A) * ...
    [b(1) 0 b(2) 0 b(3) 0;
     0 c(1) 0 c(2) 0 c(3);
     c(1) b(1) c(2) b(2) c(3) b(3)];

Ke = t * A * (B' * mat.D * B);
end
