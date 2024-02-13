%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = unifrnd(-3,3,1,1000);
s2 = unifrnd(-2,2,1,1000);
s = [s1; s2];
A = [0.6 0.8;0.8 -0.6];
x = A*s;
x1 = x(1,:);
x2 =x(2,:);

%% Question1 %%
scatter(x1,x2);
hold on

%% Question 2 %%
[aLowerSide, aLeftSide] = findMixture(x1, x2);
figure

%% Question 3 %%
x1noisy = awgn(x1,20);
x2noisy = awgn(x2,20);
scatter(x1noisy,x2noisy);
hold on
[aLowerSideNoisy, aLeftSideNoisy] = findMixture(x1noisy, x2noisy);
figure

%% Question 4 %%
histogram(x1,20);
figure

%% Question 5 %%
histogram(x2,20);
figure

%% Question 6 %%


%% PART B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1 %%
[X1,X2] = meshgrid(-10:0.1:10);
f = X1.^2 + X2.^2 - 4.*X1 - 6.*X2 + 13 + X1.*X2;
mesh(X1,X2,10*log10(f));
figure

%% Question 2 %%
contour(X1,X2,10*log10(f),'ShowText','on');
figure

%% Question 5 %% 
[minPoint1,iter1,Ys1,flag1] = steepestDescend(0.1);
[minPoint2,iter2,Ys2,flag2] = steepestDescend(0.01);

subplot(1,2,1);
plot(Ys1,'b--o');
title('Mu=0.1')

subplot(1,2,2);
plot(Ys2,'b--o');
title('Mu=0.01');

%% Question 6 %%
[minPoint3,iter3,Ys3,flag3] = newton();

%% Question 7 %%
x0 = 6;
iterMax = 25;
[points,minPointX,minPointY] = alterMinimization(iterMax,x0);
figure
plot(points(1,:),points(2,:),'b--o');
hold on
contour(X1,X2,10*log10(f),'ShowText','on');

%% Question 8 %%
Mu = 0.1;
[minPointCond,iterCond,YsCond] = steepestDescendCond(Mu);

%%