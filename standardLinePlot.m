function [] = standardLinePlot(hits, misses, tim)
    hold off
    y = mean(hits);
    plot(tim,y, ...
        'color', 'blue', 'linewidth', 2)
    stdE = std(hits) ./ sqrt(size(hits,1)); 
    y = [y+stdE flip(y)-flip(stdE)]; 
    x = [tim flip(tim)]; 
    hold on 
    fill(x, y, 'blue', 'facealpha', .2)


    y = mean(misses);
    plot(tim,y, ...
        'color', 'red', 'linewidth', 2)
    stdE = std(misses) ./ sqrt(size(misses,1)); 
    y = [y+stdE flip(y)-flip(stdE)]; 
    x = [tim flip(tim)]; 
    hold on 
    fill(x, y, 'red', 'facealpha', .2)

    xline(0, '--')







end