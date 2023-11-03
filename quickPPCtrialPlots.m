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
savFolder = [datPre 'PPC/'];
%% using chanfiles to look up the channel pairs that I want

% 
% fp = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished/';
chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;


chanFiles(~[chanFiles.targPair]) = [];

regions = unique([chanFiles.chan1Reg]); 

parfor reg1 = 1:length(regions)
    for reg2 = 1:length(regions)
%         if reg1>=reg2 && ~exist(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_231025\' 'PPC_' ...
%             regions{reg1} '_' regions{reg2} '.jpg'], 'file')
            try
        reg1Present = cellfun(@(x,y) strcmp(regions{reg1}, x) | ...
            strcmp(regions{reg1}, y), [chanFiles.chan1Reg], [chanFiles.chan2Reg]);
        reg2Present = cellfun(@(x,y) strcmp(regions{reg2}, x) | ...
            strcmp(regions{reg2}, y), [chanFiles.chan1Reg], [chanFiles.chan2Reg]);
        targPair = reg1Present & reg2Present; 
        curFiles = chanFiles(targPair);
        
        switcher = cellfun(@(x) strcmp(x, regions{reg1}), {curFiles.chan2Reg}); 

        %pull out the aggregated over channels data where each one is: 
        %trials X time X frequency X measure (PPC, TF1, TF2)
        %latency variables are trials X 
        % measure (chan1HFB, chan2HFB, RT)
        subMiss = zeros([length(curFiles)*80, 181, 10, 3]);
        smLat = zeros(length(curFiles)*80, 3); 
        smi = 1; 
        subHit = zeros([length(curFiles)*80, 181, 10, 3]);
        shLat = zeros(length(curFiles)*80, 3); 
        shi = 1; 
        retMiss = zeros([length(curFiles)*80, 161, 10, 3]);
        rmLat = zeros(length(curFiles)*80, 3); 
        rmi = 1; 
        retHit = zeros([length(curFiles)*80, 161, 10, 3]);
        rhLat = zeros(length(curFiles)*80, 3); 
        rhi = 1; 
        for ii = 1:length(curFiles)
            ii
            ch1 = split(curFiles(ii).name, '_'); 
            ch2 = split(curFiles(ii).name2, '_'); 
            subID = ch1{2}; 
            ch1 = split(ch1{3}, '.'); 
            ch2 = split(ch2{3}, '.'); 
            ch1 = str2num(ch1{1}); 
            ch2 = str2num(ch2{1}); 
            if ch1>ch2
                try
                fn = [subID '_' num2str(ch2) '_' num2str(ch1) '_' regions{reg1} '_' regions{reg2} '.mat'];
                pairDat = load([savFolder fn]).pairDat; 
                catch
                fn = [subID '_' num2str(ch2) '_' num2str(ch1) '_' regions{reg2} '_' regions{reg1} '.mat'];
                pairDat = load([savFolder fn]).pairDat; 
                end
            else

                try
                fn = [subID '_' num2str(ch1) '_' num2str(ch2) '_' regions{reg1} '_' regions{reg2} '.mat'];
                pairDat = load([savFolder fn]).pairDat; 
                catch
                fn = [subID '_' num2str(ch1) '_' num2str(ch2) '_' regions{reg2} '_' regions{reg1} '.mat'];
                pairDat = load([savFolder fn]).pairDat; 
                end

            end
            try
            encMissRT = pairDat.encInfo(pairDat.use & pairDat.misses, 4); 
            encHitRT = pairDat.encInfo(pairDat.use & pairDat.hits, 4); 
            retMissRT = pairDat.retInfo(pairDat.retInfo(:,1)==2, 3);  
            retHitRT = pairDat.retInfo(pairDat.retInfo(:,1)==1, 3); 
            if pairDat.react1 && pairDat.react2
                if switcher(ii) %reverse the latencies 
                [subMiss, smi, smLat] = updatePPCvar(subMiss, pairDat.subMissPPC, smi, ...
                    smLat, pairDat.subMissLat2, pairDat.subMissLat1, ...
                    pairDat.subMissTF2, pairDat.subMissTF1, ...
                    encMissRT); 
                [subHit, shi, shLat] = updatePPCvar(subHit, pairDat.subHitPPC, shi, ...
                    shLat, pairDat.subHitLat2, pairDat.subHitLat1, ...
                    pairDat.subHitTF2, pairDat.subHitTF1, ...
                    encHitRT); 
                [retMiss, rmi, rmLat] = updatePPCvar(retMiss, pairDat.miss_onPPC, rmi, ...
                    rmLat, pairDat.retMissLat2, pairDat.retMissLat1, ...
                    pairDat.miss_onTF2, pairDat.miss_onTF1, ...
                    retMissRT); 
                [retHit, rhi, rhLat] = updatePPCvar(retHit, pairDat.hit_onPPC, rhi, ...
                    rhLat, pairDat.retHitLat2, pairDat.retHitLat1, ...
                    pairDat.hit_onTF2, pairDat.hit_onTF1, ...
                    retHitRT);   
                else
                [subMiss, smi, smLat] = updatePPCvar(subMiss, pairDat.subMissPPC, smi, ...
                    smLat, pairDat.subMissLat1, pairDat.subMissLat2, ...
                    pairDat.subMissTF1, pairDat.subMissTF2, ...
                    encMissRT); 
                [subHit, shi, shLat] = updatePPCvar(subHit, pairDat.subHitPPC, shi, ...
                    shLat, pairDat.subHitLat1, pairDat.subHitLat2, ...
                    pairDat.subHitTF1, pairDat.subHitTF2, ...
                    encHitRT);  
                [retMiss, rmi, rmLat] = updatePPCvar(retMiss, pairDat.miss_onPPC, rmi, ...
                    rmLat, pairDat.retMissLat1, pairDat.retMissLat2, ...
                    pairDat.miss_onTF1, pairDat.miss_onTF2, ...
                    retMissRT); 
                [retHit, rhi, rhLat] = updatePPCvar(retHit, pairDat.hit_onPPC, rhi, ...
                    rhLat, pairDat.retHitLat1, pairDat.retHitLat2, ...
                    pairDat.hit_onTF1, pairDat.hit_onTF2, ...
                    retHitRT);  


                end
            end
            catch
            end
        end
        %cut off unfilled rows of the data matrices
        test = subMiss(:,1,1) == 0;
        subMiss(smi:end,:,:, :) = []; 
        subHit(shi:end,:,:, :) = []; 
        retMiss(rmi:end,:,:, :) = []; 
        retHit(rhi:end,:,:, :) = []; 
        smLat(smi:end, :) = []; 
        shLat(shi:end, :) = []; 
        rmLat(rmi:end, :) = []; 
        rhLat(rhi:end, :) = []; 
        %eliminate noise trials (participant's reaction time is WAY too
        %long)
        subMiss(smLat(:,1)==-1, :,:, :) = []; 
        smLat(smLat(:,1)==-1,:, :) = []; 
        subHit(shLat(:,1)==-1, :,:, :) = []; 
        shLat(shLat(:,1)==-1,:, :) = []; 
        retMiss(rmLat(:,1)==-1, :,:, :) = []; 
        rmLat(rmLat(:,1)==-1,:, :) = []; 
        retHit(rhLat(:,1)==-1, :,:, :) = []; 
        rhLat(rhLat(:,1)==-1,:, :) = []; 


        %figure for regional TF low frequency
        hits = subHit(:, :, :, 2); 
        misses = subMiss(:, :, :, 2); 
        hitLat = shLat(:,1); 
        missLat = smLat(:,1); 
        tim = pairDat.encTim; 
        phase = 'encode'; 
        regNam = regions{reg1}; 
        fi = [1:4]; 
        RT = shLat(:,3); 
        RT2 = smLat(:,3); 
        quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
            phase, regNam, fi, RT, RT2)

        hits = retHit(:, :, :, 2); 
        misses = retMiss(:, :, :, 2); 
        hitLat = rhLat(:,1); 
        missLat = rmLat(:,1); 
        tim = pairDat.retTim; 
        phase = 'retrieve'; 
        regNam = regions{reg1}; 
        fi = [1:4]; 
        RT = rhLat(:,3); 
        RT2 = rmLat(:,3); 
        quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
            phase, regNam, fi, RT, RT2)

        hits = subHit(:, :, :, 2); 
        misses = subMiss(:, :, :, 2); 
        hitLat = shLat(:,1); 
        missLat = smLat(:,1); 
        tim = pairDat.encTim; 
        phase = 'encode'; 
        regNam = regions{reg1}; 
        fi = [5:8]; 
        RT = shLat(:,3); 
        RT2 = smLat(:,3); 
        quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
            phase, regNam, fi, RT, RT2)

        hits = retHit(:, :, :, 2); 
        misses = retMiss(:, :, :, 2); 
        hitLat = rhLat(:,1); 
        missLat = rmLat(:,1); 
        tim = pairDat.retTim; 
        phase = 'retrieve'; 
        regNam = regions{reg1}; 
        fi = [5:8]; 
        RT = rhLat(:,3); 
        RT2 = rmLat(:,3); 
        quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
            phase, regNam, fi, RT, RT2)

%         %figure for PPC
%         figure('visible', false, 'position', [0,0,1000,1000])
%         subplot 241
%         hold off
%         histogram(shLat(:,1), [0:100:3000], 'facecolor', 'green',...
%              'normalization', 'probability')
%         hold on 
%         histogram(shLat(:,2), [0:100:3000], 'facecolor', 'k',...
%              'normalization', 'probability')
%         legend({regions{reg1}, regions{reg2}})
%         title('subsequent hit')
%         ylabel('proportion of trials')
%         xlabel('time (ms)')
% 
%         subplot 242
%         histogram(shLat(:,1) - shLat(:,2), [-2500:50:2500], 'facecolor', ...
%             'blue', 'normalization', 'probability')
%         xline(0, 'linewidth', 4, 'linestyle', '--')
%         xlabel([regions{reg1} ' leads    ' regions{reg2} ' leads'])
%         title('within trial peak dif')
%         prop = sum(shLat(:,1) - shLat(:,2) < 0) / length(shLat(:,1));
%         text(-1500, .035, num2str(round(prop, 2)) )
%         text(1000, .035, num2str(1-round(prop, 2)) )
%         ylim([0,.05])
% 
%         subplot 243
%         hold off
%         histogram(rhLat(:,1), [0:100:2500], 'facecolor', 'green',...
%              'normalization', 'probability')
%         hold on 
%         histogram(rhLat(:,2), [0:100:2500], 'facecolor', 'k',...
%              'normalization', 'probability')
%         legend({regions{reg1}, regions{reg2}})
%         title('retrieval hit')
%         ylabel('proportion of trials')
%         xlabel('time (ms)')
% 
%         subplot 244
%         histogram(rhLat(:,1) - rhLat(:,2), [-2500:50:2500], 'facecolor', ...
%             'blue', 'normalization', 'probability')
%         xline(0, 'linewidth', 4, 'linestyle', '--')
%         xlabel([regions{reg1} ' leads    ' regions{reg2} ' leads'])
%         title('within trial peak dif')
%         prop = sum(rhLat(:,1) - rhLat(:,2) < 0) / length(rhLat(:,1));
%         text(-1500, .035, num2str(round(prop, 2)) )
%         text(1000, .035, num2str(1-round(prop, 2)) )
%         ylim([0,.05])
% 
%         subplot 245
%         hold off
%         histogram(smLat(:,1), [0:100:3000], 'facecolor', 'green',...
%              'normalization', 'probability')
%         hold on 
%         histogram(smLat(:,2), [0:100:3000], 'facecolor', 'k',...
%              'normalization', 'probability')
%         legend({regions{reg1}, regions{reg2}})
%         title('subsequent miss')
%         ylabel('proportion of trials')
%         xlabel('time (ms)')
% 
%         subplot 246
%         histogram(smLat(:,1) - smLat(:,2), [-2500:50:2500], 'facecolor', ...
%             'blue', 'normalization', 'probability')
%         xline(0, 'linewidth', 4, 'linestyle', '--')
%         xlabel([regions{reg1} ' leads    ' regions{reg2} ' leads'])
%         title('within trial peak dif')
%         prop = sum(smLat(:,1) - smLat(:,2) < 0) / length(smLat(:,1));
%         text(-1500, .035, num2str(round(prop, 2)) )
%         text(1000, .035, num2str(1-round(prop, 2)) )
%         ylim([0,.05])
% 
%         subplot 247
%         hold off
%         histogram(rmLat(:,1), [0:100:2500], 'facecolor', 'green',...
%              'normalization', 'probability')
%         hold on 
%         histogram(rmLat(:,2), [0:100:2500], 'facecolor', 'k',...
%              'normalization', 'probability')
%         legend({regions{reg1}, regions{reg2}})
%         title('retrieval miss')
%         ylabel('proportion of trials')
%         xlabel('time (ms)')
% 
%         subplot 248
%         histogram(rmLat(:,1) - rmLat(:,2), [-2500:50:2500], 'facecolor', ...
%             'blue', 'normalization', 'probability')
%         xline(0, 'linewidth', 4, 'linestyle', '--')
%         xlabel([regions{reg1} ' leads    ' regions{reg2} ' leads'])
%         title('within trial peak dif')
%         prop = sum(rmLat(:,1) - rmLat(:,2) < 0) / length(rmLat(:,1));
%         text(-1500, .035, num2str(round(prop, 2)) )
%         text(1000, .035, num2str(1-round(prop, 2)) )
%         ylim([0,.05])
% 
%         set(gcf, 'color', 'w')
%         export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'withinTrialLat' ...
%             regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

% 
% 
%         frex = linspace(2,10,10); 
%         fi = [5:8];
% 
% 
%         figure('visible', true, 'position', [0,0,1000,1000])
%         maxVals = zeros(4,1); 
%         tim = pairDat.encTim; 
%         subplot 421
%         [sortedLat, order] = sort(smLat(:,1)); 
%         imagesc(tim(20:end-20), [], squeeze(mean(subMiss(order, 20:end-20, fi), 3)))
%         caxis([.25, .6])
%         hold on 
%         plot(sortedLat, [1:length(sortedLat)], 'color', 'red')
%         title(['sub miss ' regions{reg1} ' ' regions{reg2}])
% 
%         subplot 422
%         [sortedLat, order] = sort(shLat(:,1)); 
%         imagesc(tim(20:end-20), [], squeeze(mean(subHit(order, 20:end-20, fi), 3)))
%         caxis([.15, .6])
%          hold on 
%         plot(sortedLat, [1:length(sortedLat)], 'color', 'red')
%         title(['sub hit ' regions{reg1} ' ' regions{reg2}])
% 
% 
%         subplot 423
%         idx = arrayfun(@(x) find(tim>=smLat(x,1),1), [1:length(smLat(:,1))]); 
%         test = arrayfun(@(x) squeeze(mean(subMiss(x, idx(x)-20:idx(x)+20, fi),3)),...
%             1:length(idx), 'uniformoutput', false );
%         test2 = reshape(cell2mat(test), 41, []); 
%         plot(mean(test2,2), 'color', 'red', 'linewidth', 2)
%         hold on 
%         idx = arrayfun(@(x) find(tim>=shLat(x,1),1), [1:length(shLat(:,1))]); 
%         test = arrayfun(@(x) squeeze(mean(subHit(x, idx(x)-20:idx(x)+20, fi),3)),...
%             1:length(idx), 'uniformoutput', false );
%         test2 = reshape(cell2mat(test), 41, []); 
%         plot(mean(test2,2), 'color', 'blue', 'linewidth', 2)
% 
%         subplot 425
%         tim = pairDat.retTim; 
%         [sortedLat, order] = sort(rmLat(:,1)); 
%         imagesc(tim(20:end-20), [], squeeze(mean(retMiss(order, 20:end-20, fi), 3)))
% %         caxis([.15, .6])
%         hold on 
%         plot(sortedLat, [1:length(sortedLat)], 'color', 'red')
%         title(['ret miss ' regions{reg1} ' ' regions{reg2}])
% 
%         subplot 426
%         [sortedLat, order] = sort(rhLat(:,1)); 
%         imagesc(tim(20:end-20), [], squeeze(mean(retHit(order, 20:end-20, fi), 3)))
% %         caxis([.15, .6])
%          hold on 
%         plot(sortedLat, [1:length(sortedLat)], 'color', 'red')
%         title(['ret hit ' regions{reg1} ' ' regions{reg2}])
% 
% 
%         subplot 427
%         idx = arrayfun(@(x) find(tim>=rmLat(x,1),1), [1:length(rmLat(:,1))]); 
%         test = arrayfun(@(x) squeeze(mean(retMiss(x, idx(x)-20:idx(x)+20, fi),3)),...
%             1:length(idx), 'uniformoutput', false );
%         test2 = reshape(cell2mat(test), 41, []); 
%         plot(mean(test2,2), 'color', 'red', 'linewidth', 2)
%         hold on 
%         idx = arrayfun(@(x) find(tim>=rhLat(x,1),1), [1:length(rhLat(:,1))]); 
%         test = arrayfun(@(x) squeeze(mean(retHit(x, idx(x)-20:idx(x)+20, fi),3)),...
%             1:length(idx), 'uniformoutput', false );
%         test2 = reshape(cell2mat(test), 41, []); 
%         plot(mean(test2,2), 'color', 'blue', 'linewidth', 2)
       
%         set(gcf, 'color', 'w')
%         export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_231025\' 'PPC_TRIAL' ...
%             regions{reg1} '_' regions{reg2} '.jpg'], '-r300')
            catch
                disp('could not make figure')
            end
        

%         end

    end

end









