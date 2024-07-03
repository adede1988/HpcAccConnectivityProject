function [indexVals] = getScaledIndex(vals1, vals2)
   indexVals =  (vals1 - vals2) ./ (vals1 + vals2);




end