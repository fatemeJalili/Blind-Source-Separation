clc; clear all; close all;

%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('hw4-X1.mat')
load('hw4-X2.mat')

%% Question 1
subplot(2,1,1);
plot(X1(1,:));
hold on
plot(X1(2,:));
title('Question#1 - Observations');

s_est1 = MutUncoSolver(X1 , 35 , 75);

subplot(2,1,2);
plot(s_est1(1,:));
hold on
plot(s_est1(2,:));
title('Question#1 - sources');
figure 

%% Question 2
f = -50:0.2:50-0.2;
plot(f , abs(fftshift(fft(s_est1(1,:)))));
hold on 
plot(f , abs(fftshift(fft(s_est1(2,:)))));
title('Question#2 - Fourier transform of sources');
legend('s1' , 's2');
figure

%% Question 3
subplot(2,1,1);
plot(X2(1,:));
hold on
plot(X2(2,:));
title('Question#3 - Observations');

s_est2 = MutUncoSolver(X2 , 18 , 100);

subplot(2,1,2);
plot(s_est2(1,:));
hold on
plot(s_est2(2,:));
title('Question#3 - sources');
figure 

%% Question 4
f = -50:0.2:50-0.2;
plot(f , abs(fftshift(fft(s_est2(1,:)))));
hold on 
plot(f , abs(fftshift(fft(s_est2(2,:)))));
title('Question#4 - Fourier transform of sources');
legend('s1' , 's2');

%%



