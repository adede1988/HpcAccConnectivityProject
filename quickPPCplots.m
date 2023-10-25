%% quick script to get all the PPC trial data for each region pair


%% file path management

%local paths: 

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';


%HPC paths: 

% codePre = '/projects/p31578/dtf8829/';
% datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse\export_fig_repo'])

%% using chanfiles to look up the channel pairs that I want


fp = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished/';
chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;


chanFiles(~[chanFiles.targPair]) = [];

regions = unique([chanFiles.chan1Reg]); 

for reg1 = 1:length(regions)
    for reg2 = 1:length(regions)
        if reg1>=reg2 && ~exist(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_231025\' 'PPC_' ...
            regions{reg1} '_' regions{reg2} '.jpg'], 'file')
            try
        reg1Present = cellfun(@(x,y) strcmp(regions{reg1}, x) | ...
            strcmp(regions{reg1}, y), [chanFiles.chan1Reg], [chanFiles.chan2Reg]);
        reg2Present = cellfun(@(x,y) strcmp(regions{reg2}, x) | ...
            strcmp(regions{reg2}, y), [chanFiles.chan1Reg], [chanFiles.chan2Reg]);
        targPair = reg1Present & reg2Present; 
        curFiles = chanFiles(targPair);
        subMiss = zeros([length(curFiles), 181, 20]);
        subHit = zeros([length(curFiles), 181, 20]);
        retMiss = zeros([length(curFiles), 161, 20]);
        retHit = zeros([length(curFiles), 161, 20]);
        parfor ii = 1:length(curFiles)
            ii
            ch1 = split(curFiles(ii).name, '_'); 
            ch2 = split(curFiles(ii).name2, '_'); 
            ch1 = split(ch1{3}, '.'); 
            ch2 = split(ch2{3}, '.'); 
            ch1 = str2num(ch1{1}); 
            ch2 = str2num(ch2{1}); 
            chan1Dat = load([fp curFiles(ii).name]).chanDat; 
            if sum(chan1Dat.reactiveRes==1)>0
                subMiss(ii,:,:) = squeeze(chan1Dat.ISPC.subMiss(ch2, :, :, 2)); 
                subHit(ii,:,:) = squeeze(chan1Dat.ISPC.subHit(ch2, :, :, 2)); 
                retMiss(ii,:,:) = squeeze(chan1Dat.ISPC.miss_on(ch2, :, :, 2)); 
                retHit(ii,:,:) = squeeze(chan1Dat.ISPC.hit_on(ch2, :, :, 2)); 
            end
        end
        test = subMiss(:,1,1) == 0;
        subMiss(test,:,:) = []; 
        subHit(test,:,:) = []; 
        retMiss(test,:,:) = []; 
        retHit(test,:,:) = []; 
        frex = logspace(log10(2), log10(25), 20); 
       
        figure('visible', false, 'position', [0,0,1000,1000])
        maxVals = zeros(4,1); 
        tim = chan1Dat.HFB.encMulTim; 
        subplot 221
        maxVals(1) = quickPPCplotHelp(subMiss, 'sub miss', regions, reg1, reg2,tim, frex );
        subplot 222
        maxVals(2) = quickPPCplotHelp(subHit, 'sub hit', regions, reg1, reg2,tim , frex );
        tim = chan1Dat.HFB.onMulTim;
        subplot 223
        maxVals(3) = quickPPCplotHelp(retMiss, 'ret miss', regions, reg1, reg2,tim, frex );
        subplot 224
        maxVals(4) = quickPPCplotHelp(retHit, 'ret hit', regions, reg1, reg2,tim, frex );

        for pp = 1:4
            subplot(2, 2, pp)
            caxis([0, max(maxVals)])
        end
       
        set(gcf, 'color', 'w')
        export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_231025\' 'PPC_' ...
            regions{reg1} '_' regions{reg2} '.jpg'], '-r300')
            catch
                disp('could not make figure')
            end
        

        end

    end

end









