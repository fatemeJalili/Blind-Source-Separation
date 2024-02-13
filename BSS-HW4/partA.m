clc; clear all; close all;

%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1
syms x;
Cs = sin(2*pi*x) * piecewise((0<=x)&(x<1) , 0.2 , (1<=x)&(x<2) ,0.4 , (2<=x)&(x<3) , 0.6 , (3<=x)&(x<4) , -0.1 , (4<=x)&(x<5) , -0.3);
Ds = sin(4*pi*x) * piecewise((0<=x)&(x<1) , 0.1 , (1<=x)&(x<2) ,0.3 , (2<=x)&(x<3) , -0.2 , (3<=x)&(x<4) , 0.5 , (4<=x)&(x<5) , -0.3);

t = 0:0.05:5-0.05;
s1 = double(subs(Cs , t));
s2 = double(subs(Ds , t));

subplot(2,1,1);
plot(s1);
hold on
plot(s2);
title('Question#1 - Sources');

s = [s1;s2];
A = [0.8 -0.6 ; 0.6 0.8];
X = A*s;

subplot(2,1,2);
plot(X(1,:));
hold on 
plot(X(2,:));
title('Question#1 - Observations');
figure

%% Question 2
P1 = X(:,1:20);
P2 = X(:,21:40);

Rx1 = P1 * P1.';
Rx2 = P2 * P2.';

[B_tr , gamma] = eig(Rx1,Rx2);
A_est11 = normc(inv(B_tr.'));   
A_est12 = -[A_est11(:,2) , A_est11(:,1)]; % permutation ambiguity

s_est1 = A_est12 \ X;

subplot(2,1,1);
plot(s1);
hold on
plot(s2);
title('Question#2 - Sources');

subplot(2,1,2);
plot(s_est1(1,:));
hold on 
plot(s_est1(2,:));
title('Question#2 - Estimated sources');
figure

E1 = norm(s_est1 - s , 'fro').^2 / norm(s, 'fro').^2;

%% Question 3
[UReserved,DReserved] = eig(X * X.');
D = flip(flip(DReserved,1),2);
U = flip(UReserved,2);
z = sqrtm(inv(D)) * U.' * X;

P1 = z(:,1:20);
P2 = z(:,21:40);
P3 = z(:,41:60);
P4 = z(:,61:80);
P5 = z(:,81:100);
Ps = cat(3, P1, P2, P3, P4, P5);

b1 = [-0.2*10^0.5 ; 0.2*15^0.5];
b2 = [0.2*15^0.5 ; 0.2*10^0.5];

for itr = 1:100
    Rz =  pagemtimes(Ps,'none',Ps,'transpose');
   
    Rzb2 = pagemtimes(Rz,b2);
    Rzb1 = pagemtimes(Rz,b1);
    
    R1 = sum(pagemtimes(Rzb2,'none',Rzb2,'transpose'),3);
    R2 = sum(pagemtimes(Rzb1,'none',Rzb1,'transpose'),3);
    
    [U1 , D1] = eig(R1);
    [U2 , D2] = eig(R2);
    
    v = U2(:,1);
    b2 = normc((eye(2) - b1*b1.') * v);
    b1 = normc(U1(:,1));
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

E21 = norm(s_est21 - s , 'fro').^2 / norm(s , 'fro').^2;
E22 = norm(s_est22 - s , 'fro').^2 / norm(s , 'fro').^2;
[E2 , idx] = min([E21 , E22]);
if(idx == 1)
    s_est2 = s_est21;
else
    s_est2 = s_est22;
end

subplot(2,1,1);
plot(s1);
hold on
plot(s2);
title('Question#3 - Sources');

subplot(2,1,2);
plot(s_est2(1,:));
hold on 
plot(s_est2(2,:));
title('Question#3 - Estimated sources');
figure

%% Question 4
[E3, s_est3] = uncoSolver(X , s , 5 , 100);

subplot(2,1,1);
plot(s1);
hold on
plot(s2);
title('Question#4 - Sources');

subplot(2,1,2);
plot(s_est3(1,:));
hold on 
plot(s_est3(2,:));
title('Question#4 - Estimated sources');
figure

%% Question 5
[E4_k5, ~] = uncoSolver(X , s , 5 , 100);
[E4_k4, ~] = uncoSolver(X , s , 4 , 100);
[E4_k3, ~] = uncoSolver(X , s , 3 , 100);
[E4_k2, ~] = uncoSolver(X , s , 2 , 100);
E4 = [E4_k2 E4_k3 E4_k4 E4_k5];
k = [2 3 4 5];
bar(k , E4);
ylabel('E');
xlabel('k');
title('Question#5 - Error for diffrent number of windows (k)');
figure

%% Question 6
[E5_5, ~] = uncoSolver(X , s , 5 , 5);
[E5_10, ~] = uncoSolver(X , s , 5 , 10);
[E5_15, ~] = uncoSolver(X , s , 5 , 15);
[E5_20, ~] = uncoSolver(X , s , 5 , 20);
E5 = [E5_5 E5_10 E5_15 E5_20];
snr = [5 10 15 20];
bar(snr , E5);
ylabel('E');
xlabel('snr');
title('Question#6 - Error for diffrent SNRs');

%%