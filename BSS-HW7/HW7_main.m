clc; clear ; close all;
load('hw7.mat');
%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
plot(x1)
title('x1 and Shat - Q1');
hold on 
[Shat1, impulse_train1, xhat1] = singleSBD(x1);
plot(Shat1)
legend('x1','Shat')
hold off
subplot(3,1,2)
plot(impulse_train1)
title('impulse train - Q1')
subplot(3,1,3)
plot(xhat1);
title('x1hat - Q1')
% stem(Alpha,t);

%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(3,1,1)
plot(x2)
title('x2 and Shat - Q2');
hold on 
[Shat2, impulse_train2 ,xhat2] = singleSBD(x2);
plot(Shat2)
legend('x2','Shat')
hold off
subplot(3,1,2)
plot(impulse_train2)
title('impulse train - Q2')
subplot(3,1,3)
plot(xhat2)
title('x2hat - Q2')
% stem(Alpha,t);
        
%% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
k = 5;
L = 100;
T = size(x1,2);
itr = 40;
figure
subplot(3,2,1)
plot(X(1,:))
title('x1 and Shat1 - Q3')
A = normc(randi([-1,1],2,2));
for i = 1 : itr
    Bhat = pinv(A) * X;
    [Shat1, impulse_train1,xhat1] = singleSBD(X(1,:));
    [Shat2, impulse_train2,xhat2] = singleSBD(X(2,:));
    Bhat(1,:)= xhat1;
    Bhat(2,:)= xhat2;
    A = normc(X*pinv(Bhat));
end
hold on
plot(Shat1)
legend('x1','Shat1');
hold off
subplot(3,2,2)
plot(X(2,:))
hold on
title('x2 and Shat2 - Q3')
plot(Shat2)
legend('x2','Shat2');
hold off
subplot(3,2,3)
plot(impulse_train1)
title('impulse train1 - Q3')
subplot(3,2,4)
plot(impulse_train2)
title('impulse train2 - Q3')
subplot(3,2,5)
plot(xhat1)
title('x1hat - Q3')
subplot(3,2,6)
plot(xhat2)
title('x2hat - Q3')

%% Question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1 = fftshift(fft(x1));
f = -T/2 : T/2-1;
figure
plot(f/T,abs(X1)/T);
title('X1 - Q4');
% [Shat4, impulse_train4, Xhat4] = singleSBD(X1);
% hold on
% plot(Shat4)

        
    