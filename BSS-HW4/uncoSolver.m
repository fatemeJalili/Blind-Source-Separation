function [E_avg , s_est] = uncoSolver(X , s , k , snr)
E = ones(1,100);
for i = 1:100
Y = addnoise(X , snr);
[UReserved,DReserved] = eig(Y * Y.');
D = flip(flip(DReserved,1),2);
U = flip(UReserved,2);
z = sqrtm(inv(D)) * U.' * Y;

P1 = z(:,1:20);
P2 = z(:,21:40);
P3 = z(:,41:60);
P4 = z(:,61:80);
P5 = z(:,81:100);
Ps = cat(3, P1, P2, P3, P4, P5);

b1 = [-0.2*10^0.5 ; 0.2*15^0.5];
b2 = [0.2*15^0.5 ; 0.2*10^0.5];

for itr = 1:100
    Rz =  pagemtimes(Ps(:,:,1:k),'none',Ps(:,:,1:k),'transpose');
   
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

B_est11 = [b1 b2].';
B_est12 = [b1 -b2].';

B_est21 = [B_est11(:,2) , B_est11(:,1)]; % permutation ambiguity
B_est22 = [B_est12(:,2) , B_est12(:,1)]; % permutation ambiguity

s_est11 = B_est21 * z;
s_est12 = B_est22 * z;

s_est21 = [s_est11(1,:) * max(s(1,:)) / max(s_est11(1,:)) ; s_est11(2,:) * max(s(2,:)) / max(s_est11(2,:))]; % domain ambiguity
s_est22 = [s_est12(1,:) * max(s(1,:)) / max(s_est12(1,:)) ; s_est12(2,:) * max(s(2,:)) / max(s_est12(2,:))]; % domain ambiguity

E1 = norm(s_est21 - s , 'fro').^2 / norm(s , 'fro').^2;
E2 = norm(s_est22 - s , 'fro').^2 / norm(s , 'fro').^2;
E(i) = min(E1 , E2);
if(E(i) == E2)
    s_est = s_est22;
else
    s_est = s_est21;
end
i = i+1;
end
E_avg = mean(E);
end

