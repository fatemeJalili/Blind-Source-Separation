clc; clear all; close all;
load('hw6-part3.mat')
%% PART C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1
muCo = mutualCoherence(D);

%% Question 2 ,3 , 4
N0 = 2;
%%% MOD:
[Sr1,Dr1,ErrorMOD] = MOD(X,N0);
subplot(2,1,1);
plot(ErrorMOD)
title('Representation Error of OMP');
xlabel('itereation');

%%% K-SVD:
[N, ~] = size(S);
[Sr2,Dr2,ErrorKSVD] = KSVD(X, N0, N);
subplot(2,1,2);
plot(ErrorKSVD)
title('Representation Error of KSVD');
xlabel('itereation');

%% Question 5
%%% MOD:
rateMOD = recoverRt(D,Dr1);
%%% K-SVD
rateKSVD = recoverRt(D,Dr2);


