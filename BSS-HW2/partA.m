clc; clear all; close all;
%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = unifrnd(-3,3,1,1000);
s2 = unifrnd(-2,2,1,1000);
s1 = s1 - mean(s1);
s2 = s2- mean(s2);
s = [s1; s2];
A = [1 -2;2 -1;3 -2];
x = A*s;
x1 = x(1,:);
x2 =x(2,:);
x3 =x(3,:);

%% Question1 %%
scatter3(x1,x2,x3);
figure

Rx = x*x.';
[UReserved,DReserved] = eig(Rx);
D = flip(flip(DReserved,1),2);
U = flip(UReserved,2);

%% Question 2 %%
Unew = U;
Unew(:,3) = [];
C = pinv(Unew)*A;

%% Question 3 %%
Dnew = D;
Dnew(:,3)=[];
Dnew(3,:)=[];
B = sqrtm(inv(Dnew))*(Unew.');
z = B*x;
z1 = z(1,:);
z2 =z(2,:);
scatter(z1,z2);
figure

%% Question 4 %%
[Q,G,Vtranspos] = svd(x);
rank = nnz(diag(G));

%% Question 5 %%
lamdaSum = sum(diag(D));
enrgRatio = D(1,1)/lamdaSum;

u1 = U(:,1);
xRducDim = x.'*u1;
plot(xRducDim,0);

%% 




