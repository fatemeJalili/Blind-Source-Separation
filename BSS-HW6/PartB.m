clc; clear all; close all;
%% PART B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1
N1 = 10;
D1 = zeros(2,N1);
figure
for i = 1:N1
    phi = pi / N1;
    x = cos(i*phi);
    y = sin(i*phi);
    D1(:,i) = [x ; y];
    quiver(0, 0, x, y)
    hold on
end
n = 2 : 20;
figure
plot(n , cos(pi./n));
title('Question1 - mutual coherence based on N');
xlabel('N');
ylabel('Mu');

%% Question 2
N2 = 10;
D21 = zeros(3,N2);
D22 = D21;
D23 = D21;
D2 = D21;

figure
for i = 1:N2
    phi = pi / N2;
    x = cos(i*phi);
    y = sin(i*phi);
    D21(:,i) = [x ; y ; 0];
    D22(:,i) = [x ; 0 ; y];
    D23(:,i) = [0 ; x ; y];
    D2(:,i) = normc((i-1)*(i-4)*(i-7)*(i-10) * D21(:,i) + (i-2)*(i-5)*(i-8) * D22(:,i) + (i-3)*(i-6)*(i-9) * D23(:,i));
end
starts = zeros(3,3);
quiver3(zeros(1,N2), zeros(1,N2), zeros(1,N2), D2(1,:), D2(2,:), D2(3,:))
axis equal
muCo = mutualCoherence(D2);

%%