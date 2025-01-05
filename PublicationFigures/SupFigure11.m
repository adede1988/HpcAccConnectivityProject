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







%% individual peak events


    panelDat = load([datPre 'TF_HFB/'...
        'TF_hip_sub_HFB.mat']).outDat;
    ii

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals));





  % MEAN PEAK ALIGNED BROADBAND AND HFB BAND DATA

    %HITS
    BroadBandPeaksHITS = zeros([3001, length(panelDat.hitSub)]);
    HFBpeaksHITS = zeros([3001, length(panelDat.hitSub)]); 
    for tt = 1:length(panelDat.hitSub)
       
        %load in a single channel of data: 
        if tt==1
            triali = 1; 
            if panelDat.hitChi(tt)<10
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_00' num2str(panelDat.hitChi(tt)) ...
        '.mat'];
            elseif panelDat.hitChi(tt)<100
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_0' num2str(panelDat.hitChi(tt)) ...
        '.mat']; 
            else
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_' num2str(panelDat.hitChi(tt)) ...
        '.mat']; 
            end
            disp(fn)
            chanDat = load(fn).chanDat; 

            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end

        elseif (panelDat.hitChi(tt-1) ~= panelDat.hitChi(tt)) || ...
                (~strcmp(panelDat.hitSub(tt-1),panelDat.hitSub(tt)))
            triali = 1; 
            if panelDat.hitChi(tt)<10
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_00' num2str(panelDat.hitChi(tt)) ...
        '.mat'];
            elseif panelDat.hitChi(tt)<100
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_0' num2str(panelDat.hitChi(tt)) ...
        '.mat']; 
            else
                fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.hitSub{tt} '_' num2str(panelDat.hitChi(tt)) ...
        '.mat']; 
            end
            disp(fn)
            chanDat = load(fn).chanDat; 
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
        end
        

       
        %HITS
        curTrial = hfbBand(:,hitidx);
        curTrial = curTrial(:, triali); 
        L =length(curTrial); 
        curTrial = mirrorPad(curTrial); 
        phaseTrial = angle(hilbert(curTrial));
        idx = find(tim>=HFBlatHit(triali,1),1); 
        [~, troughOff] = min(phaseTrial(L+idx-10:L+idx+10));
        idx = idx + troughOff - 10; 
        HFBpeaksHITS(:, tt) = curTrial(L+idx-1500: L+idx+1500); 
        curTrial = broadBand(:, hitidx); 
        curTrial = curTrial(:, triali); 
        curTrial = mirrorPad(curTrial); 
        BroadBandPeaksHITS(:, tt) = curTrial(L+idx-1500: L+idx+1500);

        triali = triali +1; 
    end


    allRatios = []; 

    for nn = 1:1400
        %ratio of power during vs. before the HFB peak to find 'good'
        %events:
        during = mean(abs(hilbert(HFBpeaksHITS(1400:1600,nn))));
        before = mean(abs(hilbert(HFBpeaksHITS(1200:1400,nn))));
        allRatios = [allRatios during/before];

        %plotting individual events: 
        if during/before > 2 
            fig = figure('visible', false);
            subplot 211
            plot(BroadBandPeaksHITS(:,nn) - BroadBandPeaksHITS(1501,nn),...
                'color', 'k')
            hold on 
            plot(HFBpeaksHITS(:,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,3001])
            xline(1501)

            title([num2str(round(during/before, 2)) ' during/before power'])

            subplot 212
            plot(BroadBandPeaksHITS(1200:1800,nn) - BroadBandPeaksHITS(1501,nn),...
                'color', 'k')
            hold on 
            plot(HFBpeaksHITS(1200:1800,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,600])
            xline(301)
            saveas(fig, [figSavePath 'exampleTrials_supFig11/' ...
                 'highRatio_'  panelDat.reg '_'...
                 panelDat.phase '_' num2str(nn)  '.jpg'])


        end

       
    end





    %MISSES
    BroadBandPeaksMISSES = zeros([3001, length(panelDat.missSub)]);
    HFBpeaksMISSES = zeros([3001, length(panelDat.missSub)]); 

    for tt = 1:length(panelDat.missSub)
   
        %load in a single channel of data: 
        if tt==1
            triali = 1; 
            if panelDat.missChi(tt)<10
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_00' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            elseif panelDat.missChi(tt)<100
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_0' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            else
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            end
            chanDat = load(fn).chanDat; 
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
            
        elseif (panelDat.missChi(tt-1) ~= panelDat.missChi(tt)) || ...
                (~strcmp(panelDat.missSub(tt-1),panelDat.missSub(tt)))
            triali = 1; 
            if panelDat.missChi(tt)<10
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_00' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            elseif panelDat.missChi(tt)<100
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_0' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            else
fn = [datPre 'HipSingleChans/' ...
        'chanDat_' panelDat.missSub{tt} '_' num2str(panelDat.missChi(tt)) ...
        '.mat']; 
            end
            chanDat = load(fn).chanDat; 
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
        end
        
       
       
        %misses
        curTrial = hfbBand(:,missidx);
        curTrial = curTrial(:, triali); 
         L =length(curTrial); 
        curTrial = mirrorPad(curTrial); 
        phaseTrial = angle(hilbert(curTrial));
        idx = find(tim>=HFBlatMiss(triali,1),1); 
        [~, troughOff] = min(phaseTrial(L+idx-10:L+idx+10));
        idx = idx + troughOff - 10; 
        HFBpeaksMISSES(:, tt) = curTrial(L+idx-1500: L+idx+1500); 
        curTrial = broadBand(:, missidx); 
        curTrial = curTrial(:, triali); 
        curTrial = mirrorPad(curTrial); 
        BroadBandPeaksMISSES(:, tt) = curTrial(L+idx-1500: L+idx+1500);
        triali = triali +1; 
    end

    
    fig = figure('visible', false, 'position', [0,0,600,600])
    

    y = mean(BroadBandPeaksHITS(1300:1700,:),2);
    plot( y, 'color', hitCol, 'linewidth', 4)
    hold on


    y = mean(BroadBandPeaksMISSES(1300:1700,:),2);  
    plot( y, 'color', missCol, 'linewidth', 4)

    xlim([1,401])
    xticks([1:50:401])
    xticklabels([-200:50:200])
    
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;

    saveas(fig, [figSavePath  'supFig11E_broadBandHFBpeak_' panelDat.reg '_' panelDat.phase  '.jpg'])

    close(fig)














