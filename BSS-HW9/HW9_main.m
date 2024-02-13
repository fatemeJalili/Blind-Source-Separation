clc; clear ; close all;
load('hw9.mat');
%% HW9 - ICA (deflation and equivariant) 
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

%% PART A - deflation
[Xwh, whMat]  = whiten(XNoisy);

normS = norm(S,'fro');
[n,m] = size(S);

load('B5.mat')
B1 = B5;

mu = 0.1;
figure
Xwh = 100*Xwh;

for j = 1:n
    b = [];
    i = 1;
    while true
        Y = B1 * Xwh;
        bPast = B1(j,:);
        B1(j,:) = B1(j,:) - mu * ((PSI(Y(j,:), m)* Xwh.') ./ m);
        if j ~= 1
            B1(j,:) = (eye(n,n) - B1(1:j-1,:).' * B1(1:j-1,:)) * B1(j,:).';
        end
        B1(j,:) = normr(B1(j,:));
        
        b(i) = B1(j,:)*bPast.';
        i = i+1;
        
        if B1(j,:)*bPast.' > 1-1e-6
            break
        end
    end
    plot(b)
    hold on
    ylim([0 1.25])
end
title('PartA - convergance plot of 3 rows')
legend('b1*b1^T', 'b2*b2^T', 'b3*b3^T')
hold off 

Shat1 = B1*Xwh;
Shat1 = disambiguation(S,Shat1);

% for i = 1:itr
%     Y = B1 * Xwh;
%     for j = 1:n
%         B1(j,:) = B1(j,:) - mu * ((PSI(Y(j,:), m)* Xwh.') ./ m);
%         if j ~= 1
%             B1(j,:) = (eye(n,n) - B1(1:j-1,:).' * B1(1:j-1,:)) * B1(j,:).';
%         end
%         B1(j,:) = normr(B1(j,:));
%     end
%     Shat1 = B1*Xwh;
%     Shat1 = disambiguation(S,Shat1);
%     E1(i) = norm(Shat1-S, 'fro')^2/normS^2;
% end


figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat1(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part A (deflation)- real and estimated source %d' ,i))
end

permutation1 = B1*whMat*A;
Error1 = norm(Shat1-S, 'fro')^2/normS^2;

% figure 
% Error1 = E1(end);
% plot(E1)
% ylabel('Error')
% xlabel('iteration')

%% PART B - equivariant

load('B6.mat')

% B2 = rand(n,n);
B2 = B6;

itr = 500;
mu = 0.1;
E2 = ones(1,itr);
XNoisy = 100*XNoisy;

for i = 1:itr
    Y = B2 * XNoisy;
    psiY = [PSI(Y(1,:), m) ; PSI(Y(2,:), m) ; PSI(Y(3,:), m)];
    df_dt_BT = ((psiY * Y.') ./ m) - eye(n,n);
%     df_dt_BT = df_dt_BT - diag(diag(df_dt_BT));
    df_dt_BT = df_dt_BT - diag(diag(df_dt_BT)) + eye(n,n) - diag(var(Y,0,2)+ mean(Y,2));
    B2 = (eye(n,n) - mu * df_dt_BT) * B2;
    B2 = normr(B2);
    Shat2 = B2*XNoisy;
    Shat2 = disambiguation(S,Shat2);
    E2(i) = norm(Shat2-S, 'fro')^2/normS^2;
end

figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat2(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part B (equivariant)- real and estimated source %d' ,i))
end

permutation2 = B2*A;

figure 
Error2 = E2(itr);
plot(E2)
ylabel('Error')
xlabel('iteration')
title('PartB - convergance plot of Error in each iteration')

%%
