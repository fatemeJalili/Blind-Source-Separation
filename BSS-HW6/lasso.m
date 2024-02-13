function [S1] = lasso(x_noisy, D, S1, lambda)
for itr = 1 : 20
for i = 1:15
    xx = x_noisy - D * S1;
    r = xx + D(:,i) * S1(i);
    s1 = r.' * D(:,i) - lambda/2;
    s2 = r.' * D(:,i) + lambda/2;
    if (s1 > 0)
        S1(i) = s1;
    else
        if (s2 < 0)
            S1(i) = s2;
        else
            S1(i) = 0;
        end
    end
end
end
end