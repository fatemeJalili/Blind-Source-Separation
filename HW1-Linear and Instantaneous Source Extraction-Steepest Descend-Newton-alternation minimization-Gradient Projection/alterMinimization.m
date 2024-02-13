function [points,px1,px2] = alterMinimization(iterMax,x0)
iter = 0;

f = @(X1,X2)(X1^2 + X2^2 - 4*X1 - 6*X2 + 13 + X1*X2);
syms X1;
syms X2;
points = [6;6];

while (iter <= iterMax)
    px2 = double(solve (diff(f(x0,X2))== 0));
    points = [points,[x0;px2]];
    px1 = double(solve (diff(f(X1,px2))== 0));
    x0 = px1;
    points = [points,[px1;px2]];
    iter = iter + 1;
end
end