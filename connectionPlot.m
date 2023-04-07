function [] = connectionPlot(allDat, cnd1, cnd2, measure, epoch)

frex = logspace(log10(2),log10(80),100);
figure
subplot 311

test = find(cellfun(@(x) ~isempty(x), {allDat.(cnd1)}));

temp = allDat(test); 

temp = {temp.(cnd1)}; 
out1 = temp{1}; 
dimVal = length(size(out1)); 
for ii = 2:length(temp)
    out1 = cat(dimVal+1, out1, temp{ii} );
end

% out1 = sum(out1>.01, dimVal+1) ./ length(test); 
out1 = median(out1, dimVal+1, 'omitnan'); 

imagesc(imgaussfilt(out1(:,:,1)', 1))
set(gca, 'YDir', 'normal')
caxis([0,.3])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([20:20:100])
yticklabels(round(frex([20:20:100])))

title(cnd1, 'interpreter', 'none')


subplot 312


test = find(cellfun(@(x) ~isempty(x), {allDat.(cnd2)}));

temp = allDat(test); 

temp = {temp.(cnd2)}; 
out2 = temp{1}; 
dimVal = length(size(out2)); 
for ii = 2:length(temp)
    out2 = cat(dimVal+1, out2, temp{ii} );
end

% out2 = sum(out2>.01, dimVal+1) ./ length(test); 
out2 = median(out2, dimVal+1, 'omitnan'); 

imagesc(imgaussfilt(out2(:,:,1)',1))
set(gca, 'YDir', 'normal')
caxis([0,.3])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([20:20:100])
yticklabels(round(frex([20:20:100])))

title(cnd2, 'interpreter', 'none')

subplot 313



imagesc(imgaussfilt(out1(:,:,1)' - out2(:,:,1)',1))
set(gca, 'YDir', 'normal')
caxis([-.2,.2])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([20:20:100])
yticklabels(round(frex([20:20:100])))

title('difference')


end