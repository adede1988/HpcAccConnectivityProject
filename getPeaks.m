function [peaks] = getPeaks(values, lats, tim)


w = 20; 

Idx = arrayfun(@(x) find(x==tim), lats); 
Idx(Idx>length(tim)-w) = length(tim) - w; 
peaks = arrayfun(@(x) values(Idx(x)-w:...
                                 Idx(x)+w,x), [1:length(Idx)], ...
                                 'uniformoutput', false);
peaks = cell2mat(peaks);







end