clc; clear all; close all;

%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('hw3-1.mat');
TrainData_class1 = TrainData_class1 - mean(TrainData_class1,2);
TrainData_class2 = TrainData_class2 - mean(TrainData_class2,2);
TestData = TestData - mean(TestData,2);

%% Question 1
R1s = pagemtimes(TrainData_class1,'none',TrainData_class1,'transpose');
R2s = pagemtimes(TrainData_class2,'none',TrainData_class2,'transpose');

R1 = mean(R1s,3);
R2 = mean(R2s,3);

[UReserved,DReserved] = eig(R1,R2);
D = flip(flip(DReserved,1),2);
U = flip(UReserved,2);
W = normc(U);

W1 = W(:,1);
W30 = W(:,30);

E1_c1 = W1.'*TrainData_class1(:,:,49);
E1_c2 = W1.'*TrainData_class2(:,:,49);
E30_c1 = W30.'*TrainData_class1(:,:,49);
E30_c2 = W30.'*TrainData_class2(:,:,49);

var1_c1 = var(E1_c1);
var1_c2 = var(E1_c2);
var30_c1 = var(E30_c1);
var30_c2 = var(E30_c2);

subplot(1,2,1);
plot(E1_c1);
hold on 
plot(E1_c2);
legend('class#1','class#2');
title('test 49 - w1');
subplot(1,2,2);
plot(E30_c1);
hold on 
plot(E30_c2);
legend('class#1','class#2');
title('test 49 - w30');
figure

%% Question 2
subplot(1,2,1);
bar(abs(W1));
title('|W1|->class1');
subplot(1,2,2);
bar(abs(W30));
hold on 
title('|W30|->class2'); 
figure

%% Question 3
Wcsp = W;
Wcsp(:,8:23) = [];
x1 = var((pagemtimes(Wcsp, 'transpose', TrainData_class1, 'none')),0,2);
x2 = var((pagemtimes(Wcsp, 'transpose', TrainData_class2, 'none')),0,2);
x1 = reshape(x1,14,60);
x2 = reshape(x2,14,60);
Mu1 = mean(x1,2);
Mu2 = mean(x2,2);
MuDiff = Mu1-Mu2;
Sigma1 = (x1-Mu1)*(x1-Mu1).'/60;
Sigma2 = (x2-Mu2)*(x2-Mu2).'/60;
[U,D] = eig(Sigma1+Sigma2,MuDiff*MuDiff.');
Wlda = normc(U(:,1));
Mu11 = Wlda.'*Mu1;
Mu22 = Wlda.'*Mu2;
c = 0.5*(Mu11+Mu22);

%% Question 4
x = var((pagemtimes(Wcsp, 'transpose', TestData, 'none')),0,2);
x = reshape(x,14,40);
y = Wlda.'*x;
yDiff = y - c;
labelF = (yDiff > 0)+1;

%% Question 5
scatter(1:40,TestLabel,80);
hold on
scatter(1:40,labelF,30,'filled');
legend('real class','estimated class');
xlabel('test number');
ylabel('label');
ylim([0.5 2.5]);
xlim([0 41])

%%