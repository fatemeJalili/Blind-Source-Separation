function [Sr,Dr,Error] = KSVD(X,N0,N)
tic
itr = 40;
Dr = normc(randn(10, 40));
Sr = zeros(40, 1500);
Error = ones(1,itr);

for iter = 1:itr
  for i = 1:1500
      Sr(:,i) = OMP(X(:,i), Dr, N0);
  end
  for j = 1:N
      Xr = X - Dr(:, (1:N)~=j) * Sr((1:N)~=j, :);
      Xr(:, Sr(j,:) == 0) = [];
      [U,Sigma,V] = svd(Xr);
      sigMax = max(max(Sigma));
      Dr(:, j) = U(:, 1);
      VT = V.';
      Sr(j,Sr(j,:)~=0) = sigMax * VT(1,:);
  end
  Error(iter) = (norm(X - Dr * Sr,'fro')).^2;
end
toc
end