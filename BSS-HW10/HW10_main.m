clc; clear ; close all;
load('hw10.mat');
%% HW10 - ICA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% PART A - Gradient Projection (GP)
[Xwh, whMat]  = whiten(XNoisy);

Xwh = 20*Xwh;

normS = norm(S,'fro');
[n,m] = size(S);

load('B5.mat')
B1 = B5;
B2 = B1;
B3 = B1;
B4 = B1;

mu = 0.05;
figure

for j = 1:n
    b = [];
    i = 1;
    while true
        Y = B1 * Xwh;
        bPast = B1(j,:);
        kurty = mean(Y(j,:).^4) - 3*mean(Y(j,:).^2)^2;
        df_db = sign(kurty)*(mean(Xwh .* (Y(j,:).^3), 2) - 3*B1(j,:).');
        B1(j,:) = (B1(j,:).' + mu * df_db).';
        if j ~= 1
            B1(j,:) = (eye(n,n) - B1(1:j-1,:).' * B1(1:j-1,:)) * B1(j,:).';
        end
        B1(j,:) = normr(B1(j,:));
        
        b(i) = B1(j,:)*bPast.';
        i = i+1;
        
        if abs(B1(j,:)*bPast.') > 1-1e-10
       
            break
        end
    end
    plot(b)
    hold on
    ylim([0 1.25])
end
title('PartA(GP) - convergance plot of 3 rows')
legend('b1*b1^T', 'b2*b2^T', 'b3*b3^T')
hold off 

Shat1 = B1*Xwh;
Shat1 = disambiguation(S,Shat1);


figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat1(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part A (GP)- real and estimated source %d' ,i))
end

permutation1 = B1*whMat*A;
Error1 = norm(Shat1-S, 'fro')^2/normS^2;

%% PART B - Fixed Point
% B2 = rand(n,n);
figure

for j = 1:n
    b = [];
    i = 1;
    while true
        Y = B2 * Xwh;
        bPast = B2(j,:);
        B2(j,:) = mean(Xwh .* (Y(j,:).^3), 2) - 3*B2(j,:).';
        if j ~= 1
            B2(j,:) = (eye(n,n) - B2(1:j-1,:).' * B2(1:j-1,:)) * B2(j,:).';
        end
        B2(j,:) = normr(B2(j,:));
        
        b(i) = B2(j,:)*bPast.';
        i = i+1;
        
        if abs(B2(j,:)*bPast.') > 1-1e-10
            break
        end
    end
    plot(abs(b))
    hold on
    ylim([0 1.25])
end
title('Part B (FP) - convergance plot of 3 rows')
legend('b1*b1^T', 'b2*b2^T', 'b3*b3^T')
hold off 

Shat2 = B2*Xwh;
Shat2 = disambiguation(S,Shat2);


figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat2(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part B (FP)- real and estimated source %d' ,i))
end

permutation2 = B2*whMat*A;
Error2 = norm(Shat2-S, 'fro')^2/normS^2;

%% Part C - GP not sensitive to outlier
% B3 = rand(n,n);
figure

for j = 1:n
    b = [];
    i = 1;
    while true
        Y = B3 * Xwh;
        bPast = B3(j,:);
        df_db = mean(Xwh .* (Y(j,:) .* exp(-Y(j,:).^2/2)), 2);
        B3(j,:) = (B3(j,:).' + mu * df_db).';
        if j ~= 1
            B3(j,:) = (eye(n,n) - B3(1:j-1,:).' * B3(1:j-1,:)) * B3(j,:).';
        end
        B3(j,:) = normr(B3(j,:));
        
        b(i) = B3(j,:)*bPast.';
        i = i+1;
        
        if abs(B3(j,:)*bPast.') > 1-1e-10
       
            break
        end
    end
    plot(b)
    hold on
    ylim([0 1.25])
end
title('Part C (GP - not sensitive to outlier) - convergance plot of 3 rows')
legend('b1*b1^T', 'b2*b2^T', 'b3*b3^T')
hold off 

Shat3 = B3*Xwh;
Shat3 = disambiguation(S,Shat3);


figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat3(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part C (GP - not sensitive to outlier)- real and estimated source %d' ,i))
end

permutation3 = B3*whMat*A;
Error3 = norm(Shat3-S, 'fro')^2/normS^2;

%% PART D - Fixed Point  not sensitive to outlier
% B4 = rand(n,n);
figure

for j = 1:n
    b = [];
    i = 1;
    while true
        Y = B4 * Xwh;
        bPast = B4(j,:);
        B4(j,:) = mean(Xwh .* (Y(j,:) .* exp(-Y(j,:).^2/2)), 2) +...
            mean((exp(-Y(j,:).^2/2) - (Y(j,:).^2) .* (exp(-Y(j,:).^2/2))), 2) * B4(j,:).' ;
        if j ~= 1
            B4(j,:) = (eye(n,n) - B4(1:j-1,:).' * B4(1:j-1,:)) * B4(j,:).';
        end
        B4(j,:) = normr(B4(j,:));
        
        b(i) = B4(j,:)*bPast.';
        i = i+1;
        
        if abs(B4(j,:)*bPast.') > 1-1e-10
            break
        end
    end
    plot(abs(b))
    hold on
    ylim([0 1.25])
end
title('Part D (FP - not sensitive to outlier) - convergance plot of 3 rows')
legend('b1*b1^T', 'b2*b2^T', 'b3*b3^T')
hold off 

Shat4 = B4*Xwh;
Shat4 = disambiguation(S,Shat4);


figure
for i = 1 : n
    subplot(3,1,i)
    plot(Shat4(i,:))
    hold on 
    plot(S(i,:))
    legend('Shat','S')
    hold off
    title(sprintf('Part D (FP - not sensitive to outlier) - real and estimated source %d' ,i))
end

permutation4 = B4*whMat*A;
Error4 = norm(Shat4-S, 'fro')^2/normS^2;

%%


