%% Part B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1%%
k = 2*pi*150e6/3e8;
n = [0;1;2;3;4;5;6;7;8;9];
t = 0:10e-6:10e-3;
A = @(theta1, theta2) [exp(-1j*k*n*sin(theta1*pi/180)),exp(-1j*k*n*sin(theta2*pi/180))];
s = @(f1,f2) [exp(2j*pi*f1*t) ; exp(2j*pi*f2*t)];
x = A(10,20)*s(20000,10000);
xNoisy = x + wgn(10,1001,1,'linear');

%% Question 2 %%
[U,G,Vtranspos] = svd(xNoisy);
A1 = @(theta)exp(-1j*k*n*sin(theta*pi/180));
Unull = U(:,3:10);
Usig = U(:,1:2);
theta = -90:0.5:90;
f1 = vecnorm(((A1(theta)'*Usig)).');  
plot(theta,f1);
figure

%% Question 3 %%
f2 = vecnorm(((A1(theta)'*Unull)).');
plot(theta,f2);
figure

%% Question 4 %%
f = 0:100:50000;
s1 = exp(2j*pi*f.'*t);
Vnull = Vtranspos(3:1001,:);
Vsig = Vtranspos(1:2,:);
f3 = vecnorm((s1*(Vsig.')).'); 
plot(f,f3);
figure

%% Question 5 %%
f4 =  vecnorm((s1*(Vnull.')).');
plot(f,f4);
%%
