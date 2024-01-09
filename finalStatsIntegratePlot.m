%local paths: 


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
HFBdat = load([datPre 'HFB_KEY_STATS\hip.mat']).HFBdat; 
regions = {HFBdat.aggTargs.lab}; 
regions(2) = []; 




%% TF final evaluation 
%it's assumed here that the TF quest pipeline has already been run and that
%the outputs are available 

%TF power files
outStatFiles = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'all.mat'));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'phase'));
outStatFiles(test) = []; 
headFiles = dir([datPre 'TF_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, 'all.mat'));
headFiles(~test) = []; 

%TF phase files
outStatFilesPhase = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'NEWphase'));
outStatFilesPhase(~test) = []; 

%HFB files
outStatFilesHFB = dir([datPre 'HFB_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesHFB.name}, '.mat'));
outStatFilesHFB(~test) = []; 
headFilesHFB = dir([datPre 'HFB_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFilesHFB.name}, '.mat'));
headFilesHFB(~test) = []; 

%record all the significant frequencies across regions/conditions
%region X freq X sub/ret X  power/ITPC X hit/miss
HFBImageidx = cell(9, 2, 2, 2, 2);

parfor reg = 1:length(regions)
    reg
%     allOut = cell(2,2,2,2); 
%     outDat = TFfinalintegrate(reg, headFiles, outStatFiles, regions, ...
%         'sub', outStatFilesPhase);
%     allOut( :, 1, :, :) = outDat; 
%     outDat2 = TFfinalintegrate(reg, headFiles, outStatFiles, regions, ...
%         'ret', outStatFilesPhase);
%     allOut( :, 2, :, :) = outDat2; 
%     HFBImageidx(reg, :, :, :, :) = allOut; 
 
    HFBfinalintegrate(reg, headFilesHFB, outStatFilesHFB, regions, 'sub')
    HFBfinalintegrate(reg, headFilesHFB, outStatFilesHFB, regions, 'ret')




%     slice = sigFreqs(reg,:,:,:,:); 
%     slice = GEDbandDiscovery(reg, headFiles, outStatFiles, regions, ...
%         'sub', outStatFilesPhase, slice);
%     slice = GEDbandDiscovery(reg, headFiles, outStatFiles, regions, ...
%         'ret', outStatFilesPhase, slice);
%     sigFreqs(reg, :, :, :, :) = slice;
end

%% group phase plot
MTL_1 = cat(1, HFBImageidx{[1,9], :, :, 2, 1});
MTL_2 = cat(1, HFBImageidx{[1,9], :, :, 2, 2});
PFC_1 = cat(1, HFBImageidx{[2,5,7], :, :, 2, 1});
PFC_2 = cat(1, HFBImageidx{[2,5,7], :, :, 2, 2});
figure
subplot 211
hold off
histogram(MTL_1, [-1:.1:1], 'normalization', 'probability')
hold on 
histogram(PFC_1, [-1:.1:1], 'normalization', 'probability')
title('MTL v. PFC phase reset on hit trials')
xlabel(['HFB                          neutral'...
    '                          Image'])
ylim([0, .12])
subplot 212
hold off
histogram(MTL_2, [-1:.1:1], 'normalization', 'probability')
hold on 
histogram(PFC_2, [-1:.1:1], 'normalization', 'probability')
title('MTL v. PFC phase reset on miss trials')
xlabel(['HFB                          neutral'...
    '                          Image'])
ylim([0, .12])

%power plot
MTL_1 = cat(1, HFBImageidx{[1,9], :, :, 1, 1});
MTL_2 = cat(1, HFBImageidx{[1,9], :, :, 1, 2});
PFC_1 = cat(1, HFBImageidx{[2,5,7], :, :, 1, 1});
PFC_2 = cat(1, HFBImageidx{[2,5,7], :, :, 1, 2});
figure
subplot 211
hold off
histogram(MTL_1, [-1:.1:1], 'normalization', 'probability')
hold on 
histogram(PFC_1, [-1:.1:1], 'normalization', 'probability')
title('MTL v. PFC power on hit trials')
xlabel(['HFB                          neutral'...
    '                          Image'])
ylim([0, .25])
subplot 212
hold off
histogram(MTL_2, [-1:.1:1], 'normalization', 'probability')
hold on 
histogram(PFC_2, [-1:.1:1], 'normalization', 'probability')
title('MTL v. PFC power on miss trials')
xlabel(['HFB                          neutral'...
    '                          Image'])
ylim([0, .25])


% dualSignalPlot(outStatFiles, regions, outStatFilesPhase,...
%     headFiles, reg, 'sub')
