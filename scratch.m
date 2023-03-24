



load('R:\MSS\Johnson_Lab\DATA\UCI\IR84\Tasks\MemDev\BPM\data_enc_photo_bp.mat')

%% hacked hard code for electrode labels 

data.myLabels(1) = "ACC"; 
data.myLabels(2) = "ACC";
data.myLabels(3) = "OFC"; 
data.myLabels(4) = "OFC";
data.myLabels(5) = "OFC"; 
data.myLabels(6) = "OFC"; 
data.myLabels(7) = "AI"; %anterior insula
data.myLabels(8) = "BMA"; %basomedial amyg
data.myLabels(9) = "LA"; %lateral amyg
data.myLabels(10)= "LA"; 
data.myLabels(11)= "LA"; 
data.myLabels(12)= "STS"; %superior temporal sulcus
data.myLabels(13)= "STS"; 
data.myLabels(14)= "CA1"; 
data.myLabels(15)= "STS"; 
data.myLabels(16)= "MTG"; %medial temporal gyrus
data.myLabels(17)= "MTG"; 
data.myLabels(18)= "SFS"; %superior frontal sulcus 
data.myLabels(19)= "SFS";  
data.myLabels(20)= "MFG"; %medial frontal gyrus 
data.myLabels(21)= "MFG"; 
data.myLabels(22)= "MFG"; 
data.myLabels(23)= "MFG";  
data.myLabels(24)= "OFC"; %orbital frontal cortex 
data.myLabels(25)= "OFC";
data.myLabels(26)= "OFC"; 
data.myLabels(27)= "OFC"; 
data.myLabels(28)= "BLA"; %basolateral amyg 
data.myLabels(29)= "BLA"; %orbital frontal cortex 
data.myLabels(30)= "LA";
data.myLabels(31)= "LA";
data.myLabels(32)= "MTG"; 
data.myLabels(33)= "PRh"; %perirhinal 
data.myLabels(34)= "PRh"; 
data.myLabels(35)= "AFG"; %anterior fusiform gyrus 
data.myLabels(36)= "MTG";


addpath('C:\Users\dtf8829\Documents\GitHub\SheffieldAutismBiomarkers')
addpath('C:\Users\dtf8829\Documents\GitHub\HpcAccConnectivityProject')
addpath('C:\Users\dtf8829\Documents\GitHub\myFrequentUse')
addpath('C:\Users\dtf8829\Documents\GitHub\myFrequentUse\export_fig_repo')
% addpath(genpath('G:\CODE\GEDbounds_clusterImprove\eeglab2022.0\plugins\Fieldtrip-lite20220630'))

set(0,'defaultfigurewindowstyle', 'docked')

%frequency params
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);


% [nchan, ntime] = size(data.trial{1});
% ntrial = length(data.trial);
% allTrials = cell2mat(data.trial); 
% allTrials = reshape(allTrials, [nchan, ntime, ntrial]);

% locdat = struct; 
% locdat.data = allTrials; 
% locdat.srate = data.fsample;

% test = chanPower(locdat, frex, numfrex, stds); 


trialTF = getTrialTF(data.trial, frex, numfrex, stds, data.fsample);
trialHighFreq = getTrialTF(data.trial, logspace(log10(80), log10(150), 50), 50, linspace(5,8,50), data.fsample);
%test plot to demonstrate that things look normal, 
% trial, channel, time, frequency
TFplot(squeeze(trialTF(10,2,:,:))', frex, data.time{1}*1000)

%% time frequency and phase consistency plots in low frequency
%now, to make a plot of subsequently remembered and subsequently forgotten
%for each channel 
use = data.trialinfo(:,1)==1; 
hits = data.trialinfo(:,2)==1; 
misses = data.trialinfo(:,2)==2; 
pow = abs(trialTF).^2;
tic
powZ = arrayfun(@(x) myZscore(pow(:,:,:,x)), [1:size(pow,4)], 'uniformoutput', false );
toc
test = cell2mat(powZ); 
test = reshape(test, [36,4501,100,71] ); 
powZ = permute(test, [4,1,2,3]);
phase = angle(trialTF); 
timeVals = data.time{1}; 
for chan = 1:size(trialTF,2)
    baseline = squeeze(mean(abs(trialTF(use, chan, timeVals>-.45 & timeVals<-.15,:)).^2, [1,3]));
    meanPow_hit  = squeeze(mean(powZ(use & hits, chan, :,:), 1))';
    meanPow_miss  = squeeze(mean(powZ(use & misses, chan, :,:), 1))';
  

    % plot for power
    figure
    subplot 231
    imagesc(meanPow_hit)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(meanPow_hit))), max(max(abs(meanPow_hit)))])
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit power')
   
    subplot 232
    imagesc(meanPow_miss)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(meanPow_miss))), max(max(abs(meanPow_miss)))])
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('miss power')

    subplot 233
    dif = meanPow_hit - meanPow_miss; 
    imagesc(dif)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(dif))), max(max(abs(dif)))])
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title([data.label(chan)])


    % plot for ITPC
    
    
    itpc_hit = cell2mat(arrayfun(@(x) squeeze([abs(sum(exp(1i * phase(use & hits, chan,x,:) ), 1) ./ sum(use&hits) )]),...
        [1:size(meanPow_hit,2)] , 'uniformoutput', false));
    itpc_miss = cell2mat(arrayfun(@(x) squeeze([abs(sum(exp(1i * phase(use & misses, chan,x,:) ), 1) ./ sum(use&hits) )]),...
        [1:size(meanPow_hit,2)] , 'uniformoutput', false));
    
    subplot 234
    imagesc(itpc_hit)
    set(gca, 'YDir','normal')
    colorbar
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit itpc')

    subplot 235
    imagesc(itpc_miss)
    set(gca, 'YDir','normal')
    colorbar
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('miss itpc')


    subplot 236
    dif = itpc_hit - itpc_miss; 
    imagesc(dif)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(dif))), max(max(abs(dif)))])
    yticks([10:10:100])
    yticklabels(round(frex([10:10:100])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit - miss')


    set(gcf,'color','w');
    %thicken axes: 
%     axes = gca; 
%     axes.LineWidth = 4; 

%     set(gcf, 'Position',  [200, 200, 500, 400])
    export_fig(join(["H:\My Drive\Johnson\MTL_PFC_networkFigs\" 'TF_ITPC' num2str(chan) '.jpg'],''), '-r300')



end



%% time frequency and phase consistency plots in high frequency
%now, to make a plot of subsequently remembered and subsequently forgotten
%for each channel 
use = data.trialinfo(:,1)==1; 
hits = data.trialinfo(:,2)==1; 
misses = data.trialinfo(:,2)==2; 
pow = abs(trialHighFreq).^2;

tic
powZ = arrayfun(@(x) myZscore(pow(:,:,:,x)), [1:50], 'uniformoutput', false );
toc
test2 = cell2mat(powZ); 
test = reshape(test, [36,4501,50,71] ); 
powZ = permute(test, [4,1,2,3]);
% meanPow = squeeze(mean(pow(use,:,:,:), [1,3]));
% SDPow = squeeze(std(pow(use,:,:,:), [], [1,3])); 
phase = angle(trialHighFreq); 
timeVals = data.time{1}; 
for chan = 1:size(trialHighFreq,2)
    meanPow_hit  = squeeze(mean(powZ(use&hits, chan, :,:),1))';
    meanPow_miss  = squeeze(mean(powZ(use & misses, chan, :,:),1))';
%     meanPow_miss = 10*log10(meanPow_miss ./ baseline);
    frex = logspace(log10(80), log10(150), 50); 

    % plot for power
    figure
    subplot 231
    imagesc(meanPow_hit)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(meanPow_hit))), max(max(abs(meanPow_hit)))])
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit power')
   
    subplot 232
    imagesc(meanPow_miss)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(meanPow_miss))), max(max(abs(meanPow_miss)))])
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('miss power')

    subplot 233
    dif = meanPow_hit - meanPow_miss; 
    imagesc(dif)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(dif))), max(max(abs(dif)))])
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title([data.label(chan)])


    % plot for ITPC
    
    
    itpc_hit = cell2mat(arrayfun(@(x) squeeze([abs(sum(exp(1i * phase(use & hits, chan,x,:) ), 1) ./ sum(use&hits) )]),...
        [1:size(meanPow_hit,2)] , 'uniformoutput', false));
    itpc_miss = cell2mat(arrayfun(@(x) squeeze([abs(sum(exp(1i * phase(use & misses, chan,x,:) ), 1) ./ sum(use&hits) )]),...
        [1:size(meanPow_hit,2)] , 'uniformoutput', false));
    
    subplot 234
    imagesc(itpc_hit)
    set(gca, 'YDir','normal')
    colorbar
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit itpc')

    subplot 235
    imagesc(itpc_miss)
    set(gca, 'YDir','normal')
    colorbar
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('miss itpc')


    subplot 236
    dif = itpc_hit - itpc_miss; 
    imagesc(dif)
    set(gca, 'YDir','normal')
    colorbar
    caxis([-max(max(abs(dif))), max(max(abs(dif)))])
    yticks([10:10:50])
    yticklabels(round(frex([10:10:50])))
    xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
    xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
    xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)
    title('hit - miss')


    set(gcf,'color','w');
    %thicken axes: 
%     axes = gca; 
%     axes.LineWidth = 4; 

%     set(gcf, 'Position',  [200, 200, 500, 400])
    export_fig(join(["H:\My Drive\Johnson\MTL_PFC_networkFigs\" 'TF_high_ITPC' num2str(chan) '.jpg'],''), '-r300')



end

%% get time course of high frequency activity across 80-150 Hz range
frex = logspace(log10(80),log10(150),50);
highBands = [80:10:150]; 
highAmpShifts = zeros([2, size(trialHighFreq, [2,3]), length(highBands)-1]);
highBands(1) = 79; 
for fi = 1:length(highBands)-1
    
    highAmpShifts(1, :,:,fi) = mean( powZ(use&hits,:,:,frex>highBands(fi) & frex<highBands(fi+1)) , [1,4]);
    highAmpShifts(2, :,:,fi) = mean( powZ(use&misses,:,:,frex>highBands(fi) & frex<highBands(fi+1)) , [1,4]);
end

highAmpShiftHIT = squeeze(mean(highAmpShifts(1,:,:,:), 4));
highAmpShiftMISS = squeeze(mean(highAmpShifts(2,:,:,:), 4));
[labs, chanOrder] = sort(data.myLabels);
borders = find(arrayfun(@(x) ~strcmp(labs(x), labs(x+1)), 1:35));
 borders = [0, borders, 36];

figure 
subplot 121 
imagesc(highAmpShiftHIT(chanOrder, :))
ticVals = []; 
ticLabVals = [];

for b = 2:length(borders)
    %         xline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
    yline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
    ticVals = [ticVals, (borders(b-1)+1+borders(b))/2];
    ticLabVals = [ticLabVals, labs(borders(b))];
end
clim([-3,3])
yticks(ticVals)
yticklabels(ticLabVals)
xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)

title('subsequent hit high freq power')
subplot 122 
imagesc(highAmpShiftMISS(chanOrder, :))
ticVals = []; 
ticLabVals = [];

for b = 2:length(borders)
    %         xline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
    yline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
    ticVals = [ticVals, (borders(b-1)+1+borders(b))/2];
    ticLabVals = [ticLabVals, labs(borders(b))];
end
clim([-3,3])
yticks(ticVals)
yticklabels(ticLabVals)
xticks([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)])
xticklabels(round(1000*timeVals([1:round(size(meanPow_hit,2)/6):size(meanPow_hit,2)]))); 
xline(find(timeVals==0), 'color', 'red', 'linestyle', '--', 'linewidth', 3)

title('subsequent miss high freq power')
set(gcf,'color','w');

export_fig(join(["H:\My Drive\Johnson\MTL_PFC_networkFigs\" 'highFreqAmplitude' '.jpg'],''), '-r300')

%% getting ISPC values for connectivity

%now, I'd like to get all of the condition X channel X channel X time X frequency connectivity
%values
% % allISPC = zeros(2, size(trialTF,2), size(trialTF,2), size(trialTF, 3), size(trialTF,4)); 
% % 
% % % for con = 1:2
% %     for chan1 = 1:size(trialTF,2)
% %         disp(['con: ' num2str(con) ' chan: ' num2str(chan1)])
% % %         for chan2 = 1:size(trialTF,2)
% %             for tt = 1:size(trialTF,3)
% %                 for fi = 1:size(trialTF,4)
% % %                     if chan1 > chan2 %only do one half of the symetric matrix
% %                             
% %                         allISPC(1,chan1, :, tt, fi) = arrayfun(@(x) ...
% %                             abs(sum(exp(1i* bsxfun(@minus, phase(use & hits, chan1, tt, fi), phase(use & hits, x, tt, fi))))./ ...
% %                             sum(use&hits)), [1:size(trialTF,2)], 'uniformoutput', true);
% %                         allISPC(2,chan1, :, tt, fi) = arrayfun(@(x) ...
% %                             abs(sum(exp(1i* bsxfun(@minus, phase(use & misses, chan1, tt, fi), phase(use & misses, x, tt, fi))))./ ...
% %                             sum(use&misses)), [1:size(trialTF,2)], 'uniformoutput', true);
% %                        
% % 
% % %                     end
% %                 end
% %             end
% % %         end
% %     end
% % % end

% save("C:\Users\dtf8829\Documents\GitHub\HpcAccConnectivityProject\ISPCstats.mat", "allISPC", '-v7.3')
load("C:\Users\dtf8829\Documents\GitHub\HpcAccConnectivityProject\ISPCstats.mat")

%% make plots as matrices
timBreaks = [1:100:4501];
frex = logspace(log10(2),log10(80),100);
freqBreaks = logspace(log10(2),log10(80),11);
freqBreaks(1) = freqBreaks(1)-1; 
[labs, chanOrder] = sort(data.myLabels); 
 borders = find(arrayfun(@(x) ~strcmp(labs(x), labs(x+1)), 1:35));
 borders = [0, borders, 36];
for tt = 1:length(timBreaks)-1
figure
for fi = 1:10
    subplot(5,2,fi)
    curHitISPC = squeeze(mean(allISPC(1,chanOrder,chanOrder,timBreaks(tt):timBreaks(tt+1), frex>freqBreaks(fi) & frex<=freqBreaks(fi+1) ), [4,5])) ;
    imagesc(curHitISPC + curHitISPC')
    title(['hit tt: ' num2str(tt) '; ' num2str(round(freqBreaks(fi))) ':' num2str(round(freqBreaks(fi+1))) ' Hz' ] )
    clim([.05, .5])
    ticVals = []; 
    ticLabVals = []; 
    for b = 2:length(borders)
        xline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
        yline(borders(b)+.5, '--', 'linewidth', 2, 'alpha', .5)
        ticVals = [ticVals, (borders(b-1)+1+borders(b))/2];
        ticLabVals = [ticLabVals, labs(borders(b))];
    end
    xticks(ticVals)
    xticklabels(ticLabVals)
    yticks(ticVals)
    yticklabels(ticLabVals)
    axis square

end
end


%% plot on a template brain! 
elec_mni_frv_2.elec = data.elec; 
elec_mni_frv_2.label = data.myLabels; 
frex = logspace(log10(2),log10(80),100);
surface_template_l = load('R:\MSS\Johnson_Lab\smg1656\Recons\ft_templates/surface_pial_left.mat'); % from ft
surface_template_r = load('R:\MSS\Johnson_Lab\smg1656\Recons\ft_templates/surface_pial_right.mat');

freqColors = [[1:255/10:255]; zeros(10,1)'; flip([1:255/10:255]); 255*ones(10,1)'*.75]./255;
% [brewermap(10, 'blues'), ones(10,1)*.5]';
lineWidthVals = flip([1:2:20])./2;
set(0)
for tt=1:length(timBreaks)-1
    figure
    subplot 121
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    scatter3(elecLocs(:,1), elecLocs(:,2), elecLocs(:,3), 'filled', 'color', 'red')
    text(elecLocs(:,1), elecLocs(:,2), elecLocs(:,3), [data.myLabels])
    view([260,30])
    axis([-70, 70, -50, 100, -50, 50])

    subplot 122
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    scatter3(elecLocs(:,1), elecLocs(:,2), elecLocs(:,3), 'filled', 'color', 'red')
    text(elecLocs(:,1), elecLocs(:,2), elecLocs(:,3), [data.myLabels])
    view([260,30])
    axis([-70, 70, -50, 100, -50, 50])

    for fi = 1:10
        curHitISPC = squeeze(mean(allISPC(1,:,:,timBreaks(tt):timBreaks(tt+1), frex>freqBreaks(fi) & frex<=freqBreaks(fi+1) ), [4,5])) ;
        curHitISPC = curHitISPC + curHitISPC'; 
        curMissISPC = squeeze(mean(allISPC(2,:,:,timBreaks(tt):timBreaks(tt+1), frex>freqBreaks(fi) & frex<=freqBreaks(fi+1) ), [4,5])) ;
        curMissISPC = curMissISPC + curMissISPC'; 
        for chan1 = 1:size(curHitISPC,1)
            for chan2 = 1:size(curHitISPC,1)
                if chan1>chan2
                    if curHitISPC(chan1,chan2) > .3
                        if curHitISPC(chan1,chan2) >.4
                           subplot 121
                           hold on 
                           plot3(elecLocs([chan1,chan2],1), elecLocs([chan1,chan2],2),elecLocs([chan1,chan2],3),...
                                'linewidth', lineWidthVals(fi), 'color', freqColors(:,fi))
                        end
                        if curMissISPC(chan1,chan2) >.4
                           subplot 122
                           hold on 
                           plot3(elecLocs([chan1,chan2],1), elecLocs([chan1,chan2],2),elecLocs([chan1,chan2],3),...
                                'linewidth', lineWidthVals(fi), 'color', freqColors(:,fi))

                        end
                    end
                end
            end
        end
    end
    set(gcf, 'Position',  [0, 0, 2000, 1000])
    export_fig(join(["H:\My Drive\Johnson\MTL_PFC_networkFigs\" 'brainConnectionMap' num2str(tt) '.jpg'],''), '-r300')
end

%% scratch below here for testing out phase angle difference calculations
chan1 = 4; 
chan2 = 1; 
fi = 3; 
tt = 40;

mycolors = [[1:255/sum(use&hits):255]', zeros(sum(use&hits),1), flip([1:255/sum(use&hits):255])'] ./ 255; 

figure
subplot 221
vals = phase(use & hits, chan1, tt, fi);
for ii = 1:sum(use&hits)
    polarplot(repmat(vals(ii),2,1), repmat([0 1]', 1,1), 'color', mycolors(ii,:) )
    hold on 
end
subplot 222
vals = phase(use & hits, chan2, tt, fi);
for ii = 1:sum(use&hits)
    polarplot(repmat(vals(ii),2,1), repmat([0 1]', 1,1), 'color', mycolors(ii,:) )
    hold on 
end
subplot 223
histogram(phase(use & hits, chan2, tt, fi)- phase(use & hits, chan1, tt, fi))
subplot 224
histogram(angDif)
polarplot(repmat(angDif,2,1), repmat([0 1]', sum(use&hits),1), 'blue' )

trialIDX = use & hits; 
trialIDX = find(trialIDX);

cell2mat(arrayfun(@(x) myAngDif(phase(trialIDX(x), chan1, tt, fi), ...
    phase(trialIDX(x), chan2, tt, fi)), [1:length(trialIDX)], 'uniformoutput', false ))


for ii = 1:length(trialIDX)
    figure
    subplot 211
    polarplot(repmat(phase(trialIDX(ii), chan1, tt, fi),2,1), repmat([0 1]', 1,1), 'color', 'red', 'linewidth', 3)
    hold on 
    polarplot(repmat(phase(trialIDX(ii), chan2, tt, fi),2,1), repmat([0 1]', 1,1), 'color', 'blue', 'linewidth', 3)
    hold off
    title(['trial: ' num2str(ii)])

    subplot 212
    dif1 = phase(trialIDX(ii), chan1, tt, fi) -phase(trialIDX(ii), chan2, tt, fi); 
    polarplot(repmat(dif1,2,1), ...
         repmat([0 1]', 1,1), 'color', 'red', 'linewidth', 8)
    hold on 
    dif2 = mod(phase(trialIDX(ii), chan1, tt, fi) -phase(trialIDX(ii), chan2, tt, fi), 2*pi);
    polarplot(repmat(dif2,2,1), ...
         repmat([0 1]', 1,1), 'color', 'green', 'linewidth', 4)    
    dif3 = myAngDif(phase(trialIDX(ii), chan1, tt, fi), phase(trialIDX(ii), chan2, tt, fi));
    polarplot(repmat(dif3,2,1), ...
         repmat([0 1]', 1,1), 'color', 'blue', 'linewidth', 2) 
    title([num2str(round(exp(1i*dif1   ),2 ) ) ', ' num2str(round(exp(1i*dif2   ),2 ) ) ', ' num2str(round(exp(1i*dif3   ),2 ) ) ])

end


















