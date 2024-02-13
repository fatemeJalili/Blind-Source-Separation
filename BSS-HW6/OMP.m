function s = OMP(x , D , N0)
xr = x;
[~, m] = size(D);
[~, n] = size(x);
s = zeros(m,n);
for n = 1 : ceil(N0/2)
    xrr = xr;
    dott1 = xr.' *D;
    [~ , I1] = max(abs(dott1));
    s(I1) = dott1(I1);
    xr = xr - dott1(I1) * D(:,I1);
    if( (mod(N0,2) == 1) && (n == ceil(N0/2)) )
        break
    end
    dott2 = xr.' *D;
    [~ , I2] = max(abs(dott2));
    s(I2) = dott2(I2);
    s([I1 I2]) = pinv(D(:,[I1 I2]))*xrr;
    xr = xrr-D(:,[I1 I2])*s([I1 I2]);
end
