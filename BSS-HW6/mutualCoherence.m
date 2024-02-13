function muCo = mutualCoherence(A)
[~ , n] = size(A);
muCo = zeros(n);

for i = 1:n
    for j = 1:n
        if i ~= j
            muCo(i,j) = abs(A(:,i)'*A(:,j));
        end
    end
end

muCo = max(max(muCo));
end
