function [lat] = getLatency(curChan, tim, RT)

test = mean(curChan(:,tim>=-50 & tim<=RT), 1);
testTim = tim(tim>=-50 & tim <= RT); 
test = (test - min(test)); 
test = test ./ max(test); 
%test(test<.6) = 0; 
testMean = wmean(1:length(testTim), test);

lat = testTim(round(testMean));




end