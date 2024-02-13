function [aLowerSide, aLeftSide] = findMixture(x1, x2)
k = boundary(reshape(x1,[],1),reshape(x2,[],1));
%plot(x1(k),x2(k));
%hold on
%kX1Max = k(x1(k) == max(x1(k)));
kX1Min = k(x1(k) == min(x1(k)));
kX2Max = k(x2(k) == max(x2(k)));
kX2Min = k(x2(k) == min(x2(k)));

%kUpperSide = k((x2(k) >= x2(kX1Max)) & (x1(k) >= x1(kX2Max)));
kLowerSide = k((x2(k) <= x2(kX1Min)) & (x1(k) <= x1(kX2Min)));
%kRightSide = k((x1(k) >= x1(kX2Min)) & (x2(k) <= x2(kX1Max)));
kLeftSide = k((x1(k) <= x1(kX2Max)) & (x2(k) >= x2(kX1Min)));

% plot(x1(kUpperSide),x2(kUpperSide));
% hold on
plot(x1(kLowerSide),x2(kLowerSide));
hold on
% plot(x1(kRightSide),x2(kRightSide));
% hold on 
plot(x1(kLeftSide),x2(kLeftSide));

LowerSideLine = polyfit(x1(kLowerSide),x2(kLowerSide), 1);
LeftSideLine = polyfit(x1(kLeftSide),x2(kLeftSide), 1);

aLowerSide = LowerSideLine(1);
aLeftSide = LeftSideLine(1);

end