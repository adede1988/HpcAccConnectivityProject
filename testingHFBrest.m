
 fn = 'R:\MSS\Johnson_Lab\DATA\NMH\NM06\Tasks\Rest\BPR\NM06_rest.mat';
 data = load(fn).data;

 addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject')
 addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\myFrequentUse')
 addpath 'R:\MSS\Johnson_Lab\dtf8829\GitHub\myFrequentUse\export_fig_repo'
 addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
 ft_defaults; 


 dat = cat(3, data.trial{:});

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = linspace(5, 7, highnumfrex); 
tim = -600:2399;

% allpow = arrayfun(@(x) tempFunc(squeeze(dat(x,:,:)),highfrex,...
%                         1000, tim, 5), [1:2], ...
%     'uniformoutput', false);
allpow = cell(size(dat,1)); 
parfor x = 1:size(dat,1)
    tic 

    test = getChanMultiTF(squeeze(dat(x,:,:)), highfrex, 1000, tim, 5, x); 

%     [out1, out2] = tempFunc(squeeze(dat(x,:,:)),highfrex,...
%                         1000, tim, 5);
%     allpow{x} = {out1, out2}; 
%     toc
end
test = allpow; 
test(:,2:end) = []; 
test1 = cellfun(@(x) mean(x{1},3), test, 'UniformOutput',false); 
test2 = cellfun(@(x) mean(x{2},3), test, 'UniformOutput',false); 
% savedAllpow = allpow;


allpow1 = cat(3,test1{:}); 
allpow2 = cat(3,test2{:});
tim = [-600:5:2399];

%pretend that the 3000 points are -600:2399 ms relative to a trial event
%further all behavioral responses are at t=1800

allLat1 = arrayfun(@(x) gausLat(squeeze(allpow1(:,:,x)), ...
                               tim,...
                               ones(size(allpow1,2))*1800), ...
                               [1:size(allpow1,3)], 'uniformoutput', false);
allLat2 = arrayfun(@(x) gausLat(squeeze(allpow2(:,:,x)), ...
                               tim,...
                               ones(size(allpow2,2))*1800), ...
                               [1:size(allpow2,3)], 'uniformoutput', false);


for ii = 1:size(dat,1)

w = 100; 

LatIdx = arrayfun(@(x) find(x==tim), allLat1{ii}); 
LatIdx(LatIdx>length(tim)-w) = length(tim) - w; 
LatIdx(LatIdx<=w) =  w+1; 
vals1 = squeeze(allpow1(:,:,ii)); %just look at one channel 
peaks1 = arrayfun(@(x) vals1(LatIdx(x)-w:...
                             LatIdx(x)+w, x), [1:length(LatIdx)], ...
                             'uniformoutput', false);
peaks1 = cell2mat(peaks1);

LatIdx = arrayfun(@(x) find(x==tim), allLat2{ii}); 
LatIdx(LatIdx>length(tim)-w) = length(tim) - w; 
LatIdx(LatIdx<=w) =  w+1; 
vals2 = squeeze(allpow2(:,:,ii)); %just look at one channel 
peaks2 = arrayfun(@(x) vals2(LatIdx(x)-w:...
                             LatIdx(x)+w, x), [1:length(LatIdx)], ...
                             'uniformoutput', false);
peaks2 = cell2mat(peaks2);

figure 
subplot 221
[~, order] = sort(allLat1{ii}); 
imagesc(tim, [], vals1(:,order)'); 
caxis([-10,15])
subplot 223
plot([-500:5:500], mean(peaks1,2), 'color', 'blue', 'linewidth', 2)
hold on 
plot([-500:5:500], mean(peaks2,2), 'color', 'red', 'linewidth', 2)
subplot 222
[~, order] = sort(allLat2{ii}); 
imagesc(tim, [], vals2(:,order)'); 
caxis([-10,15])
subplot 224
scatter(allLat1{ii}, allLat2{ii})

end



