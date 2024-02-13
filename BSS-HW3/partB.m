clc; clear all; close all;

%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('hw3-2.mat')
data = data - mean(data,2);

X = {templateMatrix(freq(1)),templateMatrix(freq(2)),templateMatrix(freq(3)),templateMatrix(freq(4)),templateMatrix(freq(5))};

p = ones(1,5);
labelF = ones(1,15);

for i = 1:15
    Ry = data(:,:,i)*data(:,:,i).';
    for j = 1:5
        Rx = X{j}*X{j}.';
        Rxy = X{j}*data(:,:,i).';
        Ryx = data(:,:,i)*X{j}.';
        
        Sigma1 = Rx^(-0.5)*Rxy*Ry^(-1)*Ryx*Rx^(-0.5);
        [U1,D1] = eig(Sigma1);
        c = U1(:,1);
        a = Rx^(-0.5)*c;
        
        Sigma2 = Ry^(-0.5)*Ryx*Rx^(-1)*Rxy*Ry^(-0.5);
        [U2,D2] = eig(Sigma2);
        d = U2(:,1);
        b = Ry^(-0.5)*d;
        
        p(j) = abs(a.'*Rxy*b)/sqrt(a.'*Rx*a*b.'*Ry*b);
    end
    [val, idx] = max(p);
    labelF(i) = freq(idx);
end
scatter(1:15,label,80);
hold on
scatter(1:15,labelF,30,'filled');
legend('real f','estimated f');
xlabel('test number');
ylabel('f');
xlim([0 17]);
        
%%