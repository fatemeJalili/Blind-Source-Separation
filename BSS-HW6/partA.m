clc; clear all; close all;
load('hw6-part1.mat')
%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1
Spn = linprog( ones(30,1) , [] , [] , [D,-D] , x , zeros(30,1) , []);
S1 = Spn(1:15) - Spn(16:end);

%% Question 2
thr = 0.00001;
S2 = zeros(15,1);
Error = 1;
xr = x;
while (Error > thr)
    dott = xr.' *D;
    [~ , I] = max(abs(dott));
    S2(I) = dott(I);
    if (dott(I) < 1e-3)
        S2(I)= 0;
    end
    xr = xr - dott(I) * D(:,I);
    Error = norm(xr)/norm(x);
end

%% Question 3
Spn = linprog( ones(30,1) , [] , [] , [D,-D] , x_noisy , zeros(30,1) , []);
S3 = Spn(1:15) - Spn(16:end);

%% Question 4
lambda = 0.67;
S4 = lasso(x_noisy, D, S1, lambda);

