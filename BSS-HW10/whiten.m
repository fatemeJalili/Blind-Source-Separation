function [Xwh, whMat] = whiten(X)
[U,L] = eig(X*X.');
whMat = L^(-0.5) * U.' ;
Xwh = whMat * X;
end
