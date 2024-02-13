function [s_est] = MutUncoSolver(X , thau , wndl)

[UReserved,DReserved] = eig(X * X.');
D = flip(flip(DReserved,1),2);
U = flip(UReserved,2);
z = sqrtm(inv(D)) * U.' * X;

for i = 0 : floor((length(X)-wndl)/thau)
    Ps(: , : , i+1) = z(:, (1 + i*thau) : (wndl+i*thau));
end

b1 = [-0.2*10^0.5 ; 0.2*15^0.5];
b2 = [0.2*15^0.5 ; 0.2*10^0.5];

for itr = 1:100
    Rz =  pagemtimes(z(:,1:wndl),'none',Ps,'transpose');
   
    Rzb2 = pagemtimes(Rz,b2);
    Rzb1 = pagemtimes(Rz,b1);
    
    R1 = sum(pagemtimes(Rzb2,'none',Rzb2,'transpose'),3);
    R2 = sum(pagemtimes(Rzb1,'none',Rzb1,'transpose'),3);
    
    [U1 , ~] = eig(R1);
    [U2 , ] = eig(R2);
    
    v = U2(:,1);
    b1 = normc(U1(:,1));
    b2 = normc((eye(2) - b1*b1.') * v);
    itr = itr + 1;
end

B_est = [b1 b2].';
s_est = B_est * z;

end


