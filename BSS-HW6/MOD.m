function [Sr,Dr,Error] = MOD(X,N)
tic
itr = 40;
Dr = normc(randn(10, 40));
Sr = zeros(40, 1500);
Error = ones(1,itr);

for iter = 1:itr
  for i = 1:1500
      Sr(:,i) = OMP(X(:,i), Dr, N);
  end
  Dr = normc(X * pinv(Sr));
  Error(iter) = (norm(X - Dr * Sr,'fro')).^2;
end
toc
end