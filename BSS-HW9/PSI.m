function psi = PSI(y, m)
k = [ones(1,m) ; y ; y.^2 ; y.^3 ; y.^4 ; y.^5];
dk_dt = [zeros(1,m) ; ones(1,m) ; 2*y ; 3*y.^2 ; 4*y.^3 ; 5*y.^4];
E_dk_dt = mean(dk_dt, 2);
kkT = k*k.';
%kkT(kkT < 1e-3) = 0;
% E_kkT_inv = inv(kkT./m);
theta = (kkT./m)^-1 * E_dk_dt;
psi = theta.' * k;
end
