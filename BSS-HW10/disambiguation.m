function [Shat] = disambiguation(S,Shat) 
    corr = (Shat*S.');
    [m,n] = size(S);
    SL = zeros(m,n);
    for i = 1:m
        [~,indx] = max(abs(corr(i,:)));
        SL(indx,:) = Shat(i,:) * (corr(i,indx) /(Shat(i,:)*Shat(i,:).'));
    end
    Shat = SL;
end
