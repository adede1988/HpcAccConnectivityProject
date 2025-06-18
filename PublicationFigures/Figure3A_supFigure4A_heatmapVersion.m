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


%% Figure 3B and Supplemental Figure 4B ITPC spectra at image 

TF_files_image = dir([datPre 'TFphase_HFB/']); %switch to HFB files! 
TF_files_image(1:2) = []; 



for ii = 1:length(TF_files_image)
   
    panelDat = load([TF_files_image(ii).folder '/' ...
        TF_files_image(ii).name]).outDat;
 
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    %hit plot
%     figure
%     hold off
%     imagesc(squeeze(mean(panelDat.hits_image))')
%     xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
%     set(gca, 'ydir', 'normal')
%     caxis([1, 5])
% %     if ~ismember(regi, [2,5])
%     addRedOutline(panelDat.p_image, .05, 'white'); 
% %     end
%     ylim([0,50])
%     colorbar
% %     colormap(s2w2y)
%     xticks([19:20:139])
%     xticklabels(panelDat.tim([19:20:139]))
%     yticks([1:11:100])
%     yticklabels(round(panelDat.frex([1:11:100])))
%     set(gcf,'color','w');
%     box off;
%     ax=gca;ax.LineWidth=4;
% 
%     %miss plot
%     figure
%     hold off
%     imagesc(squeeze(mean(panelDat.misses_image))')
%     xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
%     set(gca, 'ydir', 'normal')
%     caxis([1, 5])
% %     if ~ismember(regi, [2,5])
%     addRedOutline(panelDat.p_image, .05, 'white'); 
% %     end
%     ylim([0,50])
%     colorbar
% %     colormap(s2w2y)
%     xticks([19:20:139])
%     xticklabels(panelDat.tim([19:20:139]))
%     yticks([1:11:100])
%     yticklabels(round(panelDat.frex([1:11:100])))
%     set(gcf,'color','w');
%     box off;
%     ax=gca;ax.LineWidth=4;


    %dif plot
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    imagesc(squeeze(mean(panelDat.hits_hfb))' - ...
        squeeze(mean(panelDat.misses_hfb))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    colorbar
    clim([-2, 2])

    addRedOutline(panelDat.p_hfb, .05, 'white'); 

    ylim([0,50])
%     colormap(s2w2y)
    xticks([9:12:41])
    xticklabels(panelDat.tim([9:12:41]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figSavePath 'Fig3B_' 'TF_phase_Dif_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')

    % make plots for the phase distribution! 
      stepSize = 8;
    %3 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figSavePath 'supFig4B_' panelDat.reg '_' ...
        panelDat.phase  '_phase_Hist.png'], '-transparent', '-r300')

    %6.5 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figSavePath 'supFig4B_' panelDat.reg '_' ...
        panelDat.phase  '_phase_Hist_high.png'], ...
        '-transparent', '-r300')

end
















