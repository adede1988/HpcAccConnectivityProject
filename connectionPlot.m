function [] = connectionPlot(allDat, cnd1, cnd2, measure, epoch)

frex = logspace(log10(2),log10(80),100);
figure('visible', false)
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

% imagesc(imgaussfilt(out1(:,:,3)', 1))
imagesc(out1(:,:,measure)')
set(gca, 'YDir', 'normal')
% caxis([-.5,.5])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([10:10:60])
yticklabels(round(frex([10:10:60])))
ylim([1, 60])

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

% imagesc(imgaussfilt(out2(:,:,3)',1))
imagesc(out2(:,:,measure)')
set(gca, 'YDir', 'normal')
% caxis([-.5,.5])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([10:10:60])
yticklabels(round(frex([10:10:60])))
ylim([1, 60])

title(cnd2, 'interpreter', 'none')

subplot 313



% imagesc(imgaussfilt(out1(:,:,3)' - out2(:,:,3)',1))
imagesc(out1(:,:,measure)' - out2(:,:,measure)')
set(gca, 'YDir', 'normal')
% caxis([-.5,.5])
colorbar
xticks([1.5:2:length(epoch)])
xticklabels(epoch([2:2:end]))
yticks([10:10:60])
yticklabels(round(frex([10:10:60])))
ylim([1, 60])

title('difference')



%save out the figure
cndBits = split(cnd1, '_');

if strfind(cndBits(1), 'sub')
export_fig(join(['R:\MSS\Johnson_Lab\dtf8829\figures\' 'conHeatMap_' 'sub' '_' cndBits(2) '_' cndBits(3)  '.jpg'],''), '-r300')
else
export_fig(join(['R:\MSS\Johnson_Lab\dtf8829\figures\' 'conHeatMap_' cndBits(2) '_' cndBits(3) '_' cndBits(4)  '.jpg'],''), '-r300')
end





end