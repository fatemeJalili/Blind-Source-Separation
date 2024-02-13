function [Wcsp, Wlda, c, dirc] = train2class(TrainData_class1, TrainData_class2 , accuracy_level)

TrainData_class1 = TrainData_class1 - mean(TrainData_class1,2);
TrainData_class2 = TrainData_class2 - mean(TrainData_class2,2);

R1s = pagemtimes(TrainData_class1,'none',TrainData_class1,'transpose');
R2s = pagemtimes(TrainData_class2,'none',TrainData_class2,'transpose');

R1 = mean(R1s,3);
R2 = mean(R2s,3);

[UReserved,DReserved] = eig(R1,R2);
W = normc(UReserved);
[~,indx] = sort(diag(DReserved));
%D = Dreserved(ind,ind);
W = W(:,indx);

%accuracy_level = 8;

Wcsp = W;
Wcsp(:,accuracy_level+1 : size(W, 1)-accuracy_level) = [];
absWcsp = sum(abs(Wcsp),2);
figure
bar(absWcsp)

x1 = var((pagemtimes(Wcsp, 'transpose', TrainData_class1, 'none')),0,2);
x2 = var((pagemtimes(Wcsp, 'transpose', TrainData_class2, 'none')),0,2);
x1 = reshape(x1, size(x1,1), size(x1,3));
x2 = reshape(x2, size(x2,1), size(x2,3));

figure
scatter(x1(1,:),x1(2,:))
hold on
scatter(x2(1,:),x2(2,:))

Mu1 = mean(x1,2);
Mu2 = mean(x2,2);
MuDiff = Mu1-Mu2;

Sigma1 = (x1-Mu1)*(x1-Mu1).'/size(x1,2);
Sigma2 = (x2-Mu2)*(x2-Mu2).'/size(x2,2);

[Ureserved, Dreserved] = eig(Sigma1 + Sigma2 , MuDiff * MuDiff.');
U = normc(Ureserved);
[~,indx] = sort(diag(Dreserved));
%D = Dreserved(ind,ind);
U = U(:,indx);
Wlda = U(:,end);

Mu11 = Wlda.' * Mu1;
Mu22 = Wlda.' * Mu2;
if Mu11 > Mu22
    dirc = 1;
else
    dirc = 0;
end
c = 0.5 * (Mu11 + Mu22);
end


% %% Question 4
% x = var((pagemtimes(Wcsp, 'transpose', TestData, 'none')),0,2);
% x = reshape(x,14,40);
% y = Wlda.'*x;
% yDiff = y - c;
% labelF = (yDiff > 0)+1;
% end

%%