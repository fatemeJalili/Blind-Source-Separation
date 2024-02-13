function rate = recoverRt(D,Dr) 
    success = 0;
    for i = 1:40
        [maxCor,indx] = max(abs(Dr(:,i).'*D));
        if maxCor > 0.98
            D(:,indx) = [];
            success = success + 1;
        end
    end
    rate = success / 40;
end
