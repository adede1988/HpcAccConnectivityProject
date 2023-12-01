function [] = tValplot(tvals, tim, p)


hold off
plot(tim, tvals,...
    'color', 'k', 'linewidth', 2)

tim1 = tim(p<.1); 
tim2 = tim(p<.05); 
hold on 

breaks = find(diff(tim1)>25); 
breaks = [0, breaks, length(tim1)];
for ii=1:length(breaks)-1
    plot(tim1(breaks(ii)+1:breaks(ii+1)),...
        zeros(breaks(ii+1)-breaks(ii), 1),...
    'color', 'green', 'linewidth', 2, 'linestyle', '-');
end
   
breaks = find(diff(tim2)>25); 
breaks = [0, breaks, length(tim2)];
for ii=1:length(breaks)-1
    plot(tim1(breaks(ii)+1:breaks(ii+1)),...
        zeros(breaks(ii+1)-breaks(ii), 1),...
    'color', 'red', 'linewidth', 4, 'linestyle', '-');
end




end