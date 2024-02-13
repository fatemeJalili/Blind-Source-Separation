clc; clear all; close all;
load('hw5.mat')
N0 = 3;
%% Sparce Signal Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q1 - Subset Selection 
tic;
A = cell(1,60);
for i= 1:60
    A{1,i} = D(:,i);
end
B = nchoosek(A,3);
indx = nchoosek(1:60,3);

E = ones(34220,1);
S = cell(34220,1);
for i = 1:34220
    Dp = [B{i,1},B{i,2},B{i,3}];
    S{i,1} = pinv(Dp) * x;
    E(i,1) = norm(x - Dp * S{i,1});
end
[~,I] = min(E);

non0 = indx(I,:);
ss = cell2mat(S(I));
S1 = zeros(60,1);
for i = 1:3
  S1(non0(i)) = ss(i);  
end
toc;

%% Q2 - minimize norm2
tic;
S2 = pinv(D)*x;
toc;

%% Q3 - Matching Pursuit
%%% known N0:
tic;
xr = x;
S31 = zeros(60,1);
for n = 1:N0
    dott = xr.' *D;
    [~ , I] = max(abs(dott));
    S31(I) = dott(I);
    xr = xr - dott(I) * D(:,I);
end
toc;

tic;
thr = 0.01;
S32 = zeros(60,1);
Error = 1;
xr = x;
while (Error > thr)
    dott = xr.' *D;
    [~ , I] = max(abs(dott));
    S32(I) = dott(I);
    xr = xr - dott(I) * D(:,I);
    Error = norm(xr)/norm(x);
end
toc;
% when I run this part in command window or run this section separately it gives me diffrent Elapsed time!
% Elapsed time in report is based on which command window gives
%% Q4 - Orthogonal Matching Pursuit
%%% known N0:
tic;
xr = x;
S41 = zeros(60,1);
for n = 1 : ceil(N0/2)
    xrr = xr;
    dott1 = xr.' *D;
    [~ , I1] = max(abs(dott1));
    S41(I1) = dott1(I1);
    xr = xr - dott1(I1) * D(:,I1);
    if( (mod(N0,2) == 1) && (n == ceil(N0/2)) )
        break
    end
    dott2 = xr.' *D;
    [~ , I2] = max(abs(dott2));
    S41(I2) = dott2(I2);
    S41([I1 I2]) = pinv(D(:,[I1 I2]))*xrr;
    xr = xrr-D(:,[I1 I2])*S41([I1 I2]);
end
toc;

%%% unknown N0
tic;
xr = x;
S42 = zeros(60,1);
Error = 1;
while (Error > thr)
    xrr = xr;
    dott1 = xr.' *D;
    [~ , I1] = max(abs(dott1));
    S42(I1) = dott1(I1);
    xr = xr - dott1(I1) * D(:,I1);
    dott2 = xr.' *D;
    [~ , I2] = max(abs(dott2));
    S42(I2) = dott2(I2);
    S42([I1 I2]) = pinv(D(:,[I1 I2]))*xrr;
    xr = xrr-D(:,[I1 I2])*S42([I1 I2]);
    Error = norm(xr)/norm(x);
end
toc;

%% Q5 - Basis Pursuit

tic;
Spn = linprog( ones(120,1) , [] , [] , [D,-D] , x , zeros(120,1) , []);
S5 = Spn(1:60)-Spn(61:end);
toc;
% when I run this part in command window or run this section separately it gives me 10 times lower Elapsed time!
% Elapsed time in report is based on which command window gives

%% Q6 - Iteratively Reweighted Least Square
tic;
S6 = zeros(60,1);
w = rand(1,60);
for i = 1:70
    P = D * diag(w.^-0.5);
    y = pinv(P) * x;
    S6 = diag(w.^-0.5) * y;
    for j = 1:60
        if(abs(S6(j)) > 1e-7)
            w(j) = (abs(S6(j))).^-1;
        else
            w(j) = 10e7;
        end
    end
end
S6(abs(S6) < 1e-7) = 0;
toc;

%%