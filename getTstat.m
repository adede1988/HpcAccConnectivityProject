function [tVal] = getTstat(x,y)
    [~,~,~,stats] = ttest2(x,y); 
    tVal = stats.tstat;
end