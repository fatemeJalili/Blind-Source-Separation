clc; clear ; close all;
load('hw8.mat');
%% HW8 - ICA (using argmin D_kl and MSE) %%%%%%%
X = A*S;
XNoisy = X + Noise;

subplot(3,1,1)
plot(S(1,:))
hold on 
plot(S(2,:))
plot(S(3,:))
legend('S1','S2','S3');
title('Sources')
hold off

subplot(3,1,2)
plot(X(1,:))
hold on
plot(X(2,:))
plot(X(3,:))
legend('X1','X2','X3');
title('Observations without Noise')
hold off

subplot(3,1,3)
plot(XNoisy(1,:))
hold on 
plot(XNoisy(2,:))
plot(XNoisy(3,:))
legend('XNoisy1','XNoisy2','XNoisy3');
title('Noisy Obserations')

%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normS = norm(S,'fro');
[n,m] = size(S);
y = (2*rand(1, m)-1)/10;
psi = PSI(y, m);
%B = 4*rand(n,n) -1 ;
B = [2 1 0; -1 -1 1 ; 2 0 -1];
itr = 5000;
mu = 0.01;
E = ones(1,itr);
for i = 1:itr
    Y = B * XNoisy;
    psiY = [PSI(Y(1,:), m) ; PSI(Y(2,:), m) ; PSI(Y(3,:), m)];
    df_dt = ((psiY * XNoisy.') ./ m) - inv(B.');
    B = B - mu * df_dt;
    Shat = B*XNoisy;
    E(i) = norm(Shat-S, 'fro')^2/normS^2;
end
figure
plot(Shat(1,:))
hold on 
plot(Shat(2,:))
plot(Shat(3,:))
legend('S1','S2','S3');
title('Estimated Sources')
hold off

figure
plot(S(1,:))
hold on 
plot(S(2,:))
plot(S(3,:))
plot(Shat(1,:))
plot(Shat(2,:))
plot(Shat(3,:))
legend('S1','S2','S3','Shat1','Shat2','Shat3');
title('Estimated Sources and real Sources')
hold off

figure 
Error = E(itr);
plot(E)
ylabel('Error')
xlabel('iteration')






