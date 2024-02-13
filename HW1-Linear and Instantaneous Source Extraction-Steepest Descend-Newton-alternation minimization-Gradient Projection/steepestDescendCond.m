function [x0,iter,Ys] = steepestDescendCond(Mu)
x0 = [6;6];
iter_max = 10000;
tolerance = 10e-7;
iter = 0;
grad = [1; 1]; 
xNewCond = [1; 1];

g1 = @(X1,X2)(2 * X1 - 4 + X2);
g2 = @(X1,X2)(2 * X2 - 6 + X1);
Ys = x0(2);
syms xx;
while (norm(grad,2) >= tolerance)
    grad(1,1) = g1(x0(1), x0(2));
    grad(2,1) = g2(x0(1), x0(2));
    xNew = x0 - Mu * grad; 
    xNewCond = xNew / norm(xNew,2);
    x0 = xNewCond;
    Ys = [Ys,x0(2)];
    iter = iter + 1;
    if iter > iter_max
        break
    end
end
end