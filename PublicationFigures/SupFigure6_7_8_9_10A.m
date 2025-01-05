%% MTL - PFC publication figures 

%description of the general data set 


%% set environment

%figures should not be docked: 
set(0, 'defaultfigurewindowstyle', 'normal')
%local paths: 

%edit here for local file paths: 
codePre = 'G:\My Drive\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\Pubdat\';
figDat = datPre;
figSavePath = 'R:\MSS\Johnson_Lab\dtf8829\Pubdat\FiguresOut\';

% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
% addpath([codePre 'myFrequentUse'])
% addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath([codePre 'fieldtrip-20230118'])
% ft_defaults;
regions = {'acc', 'dlPFC', 'hip', 'mtl',  'pPFC'}; 




%% set graphical parameters: 

hitCol = [88, 113, 206] ./ 256; 
missCol = [240, 172, 99] ./ 256; 
sigCol = [.5,.5,.5]; 
sigAlpha = .3; 
errAlpha = .2; 

regColors = [[204, 153, 204]', ...%ACC: dusky pink
             [145, 162, 80]', ... %dlPFC: light forest green
             [61, 187, 218]', ... %hip light blue
             [214, 77, 97]',...   %parahip: dark pink
             [51, 102, 153]']' ./ 255; %polarPFC: dark blue
            

keyRegIdx = [1,2,3,4,5]; 

b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], ...
    [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)';...
    linspace(255,0,128)']]/255;
b2w2r(129:end, 1) = 1; 
b2w2r(1:128, 3) = 1; 

linWid = 5; 

%color anchors for custom color map
s = [12, 61, 74]; 
m = [171,189,154]; 
e = [255, 255, 46]; 

s2w2y = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
         [linspace(s(2),m(2),128)'; linspace(m(2),e(2),128)'], ...
         [linspace(s(3),m(3),128)'; linspace(m(3),e(3),128)'], ...
         ] / 255;

%color anchors for custom color map
s = [85,25,86];
m = [186,21,77]; 
e = [249,205,15]; 

purpleYellow = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
         [linspace(s(2),m(2),128)'; linspace(m(2),e(2),128)'], ...
         [linspace(s(3),m(3),128)'; linspace(m(3),e(3),128)'], ...
         ] / 255;


phaseVals = {'sub', 'ret'}; 

figure
hold on 
for ii = 1:5
    scatter(ii, 2, 100, regColors(ii,:), 'filled')
    text(ii, 2.5, regions{ii})

end



%% Supplemental Figures 6 through 9 and 10A grab all connections: 

fig4Dat = dir([datPre 'connectionDat/']);  
fig4Dat(1:2) = []; 

%all connections in a reg X reg X time X hit/miss/t/p X enc/ret
allConnections = zeros(9,9,139, 20, 4, 2); %image locked
allConnections2 = zeros(9,9,41,20,4,2); %HFB locked
for ii = 1:length(fig4Dat)
    curDat = load([fig4Dat(ii).folder '/' fig4Dat(ii).name]).outDat; 
    disp(ii)
    
    
    subIDs = curDat.subVals; 
    uniIDs = unique(subIDs);
    varCodes = {'enc_image', 'ret_image', 'enc_HFB', 'ret_HFB'};
    direction = {'pos', 'neg'}; 
    for d = 1:2
    for cc = 1:4
    %do encoding image locked
    if isfield(curDat.([varCodes{cc} '_clust']), [direction{d} '_clusters'])
        clust = curDat.([varCodes{cc} '_clust']).([direction{d} '_clusters']); 
        clust(:,[clust.p]>.05) = [];
        for clu = 1:length(clust)
            if ~isempty(clust(clu).p)
            idx = clust(clu).inds; 
            L = sum(idx, 'all');
            n = size(curDat.([varCodes{cc} '_hitVals']),3);

            chanMeans = []; 
            chanMeansHit = []; 
            chanMeansMiss = []; 
            for chan = 1:n
                sliceHit = squeeze(curDat.([varCodes{cc} '_hitVals'])(:,:,chan));
                sliceMiss= squeeze(curDat.([varCodes{cc} '_missVals'])(:,:,chan));
                chanMeansHit = [chanMeansHit, mean(sliceHit(idx))];
                chanMeansMiss = [chanMeansMiss, mean(sliceMiss(idx))]; 
                slice = sliceHit - sliceMiss; 
                chanMeans = [chanMeans, mean(slice(idx))]; 
            end




% make linked boxplot 
    figure('visible', false, 'position', [0,0,600,400])
    hmSort = [zeros(1,n), ones(1,n)]; 
    b1 = boxchart(zeros(1,n), chanMeansMiss, 'boxfacecolor', missCol);
    hold on 
    b2 = boxchart(ones(1,n), chanMeansHit, 'boxfacecolor', hitCol);
    allChanMeans = [chanMeansMiss, chanMeansHit]; 
    b1.MarkerStyle = 'none'; 
    b2.MarkerStyle = 'none'; 
    b1.LineWidth = 3;  % This changes the outline of the box, but not the whiskers
    b2.LineWidth = 3;  % Same here
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
%             meanX = round(mean(Xidx(tmp), 'all')); 
%             meanY = round(mean(Yidx(tmp), 'all')); 
%             LLvals = -150:150; 
    %DRAW CUSTOM COLORED ERROR BARS: 
    IQR = prctile(chanMeansMiss, [25,75]);
    Q25 = IQR(1); 
    Q75 = IQR(2); 
    IQR = diff(IQR); 
    tmp = chanMeansMiss(chanMeansMiss>(Q25-1.5*IQR) & chanMeansMiss<(Q75+1.5*IQR));
    % Overlay custom error bars
    % For Miss group
    line([-.15,.15], [max(tmp), max(tmp)],  'LineWidth', 3, 'Color', missCol);
    line([0,0], [Q75, max(tmp)],  'LineWidth', 3, 'Color', missCol);
    
    line([-.15,.15], [min(tmp), min(tmp)],  'LineWidth', 3, 'Color', missCol);
    line([0,0], [Q25, min(tmp)],  'LineWidth', 3, 'Color', missCol);
    
    IQR = prctile(chanMeansHit, [25,75]);
    Q25 = IQR(1); 
    Q75 = IQR(2); 
    IQR = diff(IQR); 
    tmp = chanMeansHit(chanMeansHit>(Q25-1.5*IQR) & chanMeansHit<(Q75+1.5*IQR));
    % Overlay custom error bars
    % For hit group
    line([1-.15,1.15], [max(tmp), max(tmp)],  'LineWidth', 3, 'Color', hitCol);
    line([1,1], [Q75, max(tmp)],  'LineWidth', 3, 'Color', hitCol);
    
    line([1-.15,1.15], [min(tmp), min(tmp)],  'LineWidth', 3, 'Color', hitCol);
    line([1,1], [Q25, min(tmp)],  'LineWidth', 3, 'Color', hitCol);
    
    PLH = allChanMeans(hmSort==1); 
    PLM = allChanMeans(hmSort==0); 
    randVals = (rand(length(PLH),1)-.5)*.5;
    hold on 
    scatter(randVals, PLM, 10,  'blue')
    scatter(randVals+1, PLH, 10, 'blue')
    
    for ppi = 1:length(PLH)
        plot([0+randVals(ppi),1+randVals(ppi)], [PLM(ppi),PLH(ppi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
        

    end

    %make subject means
    subHits = zeros(length(uniIDs),1); 
    subMisses = zeros(length(uniIDs),1); 
    uniqueSubs = uniIDs; 
    regSubs = subIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    xlim([-.5, 1.5])
%     ylim([-.05, min([.3; max([PLH, PLM])+.05])])



set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
export_fig([figSavePath 'supFig10A_' 'ppc_' direction{d} '_' varCodes{cc} '_' ...
                curDat.reg1 '_' curDat.reg2 '_' num2str(clu) '.jpg'], '-r300')


            end
        end
    end
    end
    end




    tim = curDat.enc_image_tim; 
    reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
    reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions)); 
    %ENCODING: 
    %hit vals
    tmpMat = curDat.enc_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    %hit image figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_image_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig6A_' 'enc_' 'connectionHeatmap_image_hit_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.enc_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    %hit HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_HFB_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    pMat = curDat.enc_HFB_p;
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig7A_' 'enc_' 'connectionHeatmap_HFB_hit_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    %miss vals
    tmpMat = curDat.enc_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3); 


    %miss image figure
     %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_image_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig6B_' 'enc_' 'connectionHeatmap_image_miss_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.enc_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3);

    %miss HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_HFB_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000);  
    pMat = curDat.enc_HFB_p;
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig7B_' 'enc_' 'connectionHeatmap_HFB_miss_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    %t vals
    tmpMat = curDat.enc_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 1) = tmpMat; 

    %image t values
    %ss: 1 = HFB, 0 = image
    inMat = curDat.enc_image_tVal;  
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig6C_' 'enc_' 'connectionHeatmap_image_tVal_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')



    tmpMat = curDat.enc_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 1) = tmpMat;

    %HFB t values 
     %ss: 1 = HFB, 0 = image
    inMat = curDat.enc_HFB_tVal;  
    pMat = curDat.enc_HFB_p; 
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig7C_' 'enc_' 'connectionHeatmap_HFB_tVal_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %p vals
    tmpMat = curDat.enc_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 1) = tmpMat; 

    tmpMat = curDat.enc_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 1) = tmpMat;


    %RETRIEVAL: 
    tim = curDat.ret_image_tim; 
    %hit vals
    tmpMat = curDat.ret_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3); 

    %hit image figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_image_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig8A_' 'ret_' 'connectionHeatmap_ret_image_hit_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.ret_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3);

    %hit HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_HFB_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    pMat = curDat.ret_HFB_p;
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig9A_' 'ret_' 'connectionHeatmap_ret_HFB_hit_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %miss vals
    tmpMat = curDat.ret_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 


    %miss image figure
     %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_image_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig8B_' 'ret_' 'connectionHeatmap_ret_image_miss_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    tmpMat = curDat.ret_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 

    %miss HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_HFB_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000);  
    pMat = curDat.ret_HFB_p;
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig9B_' 'ret_' 'connectionHeatmap_ret_HFB_miss_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %t vals
    tmpMat = curDat.ret_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 2) = tmpMat; 

     %image t values
    %ss: 1 = HFB, 0 = image
    inMat = curDat.ret_image_tVal;  
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig8C_' 'ret_' 'connectionHeatmap_ret_image_tVal_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    tmpMat = curDat.ret_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 2) = tmpMat; 

     %HFB t values 
     %ss: 1 = HFB, 0 = image
    inMat = curDat.ret_HFB_tVal;  
    pMat = curDat.ret_HFB_p; 
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400]);
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    export_fig([figSavePath 'supFig9C_' 'ret_' 'connectionHeatmap_ret_HFB_tVal_'...
        regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %p vals
    tmpMat = curDat.ret_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 2) = tmpMat; 

    tmpMat = curDat.ret_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 2) = tmpMat; 

end

save([datPre 'HFBConnections.mat'], 'allConnections2')
save([datPre 'imageConnections.mat'], 'allConnections')













