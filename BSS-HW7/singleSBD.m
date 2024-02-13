function [Shat, impulse_train ,xhat] = singleSBD(x1)
k = 5;
L = 100;
T = size(x1,2);
itr = 30;
Alpha = 0.2 * ones(k,1);
Y = ones(k,L);
t = 1 : T/k : T;
for i = 1:itr
    for l = 1:k
        Y(l,:) = x1(t(l):t(l)+L-1);
    end
    Shat = (Y.' * Alpha) / (Alpha.'* Alpha) ;
    Shat = Shat/norm(Shat);
%     plot(Shat)
    [corr , lag] = xcorr(x1 , Shat);
    corr = corr(lag>=0 & lag<=T-L);
%     corr = corr(T:2*T-1);
%     hold on
%     plot(corr)
    for j = 1: k
        [Alpha(j) , t(j)] = max(corr);
        if (t(j)-L +1 < 1)
            corr(1 : t(j)+L-1) = 0;
        else
            if (t(j)+L-1 > T-L)
                corr (t(j)-L+1 : T-L) = 0;
            else
                corr (t(j)-L+1 : t(j)+L -1) = 0;
            end
        end
%         hold on
%         plot(corr)
    end
end 
xhat = zeros(1,T);
impulse_train = zeros(1,L);
for i = 1:k
    impulse_train(t(i)+L/2) = Alpha(i);
    xhat(t(i):t(i)+L-1) = Shat*Alpha(i);
end
end
