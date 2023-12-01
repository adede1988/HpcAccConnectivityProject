function [pow, pow2] = tempFunc(dat, highfrex, sr, tim, window)

    %dat needs to be in a time X trials matrix 
    [pow, mulTim, mulFrex] = ...
                getChanMultiTF(dat, highfrex, sr, tim, window); 

    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), ...
        [1, length(mulTim)] ), ...
        1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize


%% fieldtrip based 
    [spectrum, ntaper, freqoi, timeoi] = ft_specest_mtmconvol(dat', ...
        tim/1000, 'freqoi', linspace(80,150,25), ...
        'timeoi', tim(1:5:end)/1000,...
        'timwin', 5./linspace(80,150,25), ...
        'taper', 'dpss', ...
        'tapsmofrq', 0.4 *linspace(80,150,25),...
        'pad', 30);
    pow2 = abs(spectrum).^2;
    pow2 = squeeze(mean(pow2, 1, 'omitnan')); 
    pow2 = permute(pow2, [3,1,2]);
   
    pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), ...
        [1, length(timeoi)] ), ...
        1:size(pow2,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(freqoi); 
    pow2 = cell2mat(pow2); %organize
    pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []);

%     %% debug work: 
%     pow = squeeze(mean(pow,3)); 
%     pow2 = squeeze(mean(pow2,3)); 
% 
% 
%     lat = gausLat(pow, timeoi*1000, ones(size(pow,2))*1800)  ;
%     lat2 = gausLat(pow2, timeoi*1000, ones(size(pow,2))*1800)  ;
%     
%     figure
%     [~, order] = sort(lat); 
%     subplot 211
%     imagesc(pow(:,order)')
%     caxis([-10, 20])
%     [~, order] = sort(lat2); 
%     subplot 212
%     imagesc(pow2(:,order)')
%     caxis([-10, 20])
% 




end