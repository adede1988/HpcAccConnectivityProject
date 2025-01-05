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
addpath([codePre 'myFrequentUse/export_fig_repo'])
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

%% Figure 1C and 2A: methods figure, what is HFB reactive? 

test = load([datPre 'HFB_singleTrial/HFB_singleTrialmtl_sub_image.mat']).outDat;

%load in two example channels of raw data from one subject to plot
%examples of reactive and non reactive channels
reactChan = load([datPre 'chanDat_IR84_034.mat']).chanDat; 
testSub = test.hitChi(cellfun(@(x) strcmp(x, 'IR84'), test.hitSub));
nonReact = load([datPre 'chanDat_IR84_035.mat']).chanDat;

figure('position', [0,0,600,1200])
allTrials = [reactChan.HFB.subHit, reactChan.HFB.subMiss]; 
imagesc(reactChan.HFB.encMulTim, [], allTrials')
colormap(s2w2y)
clim([0,10])
colorbar
hold on 
xline(0, '--', 'linewidth', linWid, 'color', 'white')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
yticks([])
xlim([-450, 3000])
export_fig([figSavePath 'Fig1ci_reactiveMethod_yes_heat' '.jpg'], '-r300')

figure('position', [0,0,600,400])
plot(reactChan.HFB.encMulTim, mean(allTrials,2), 'linewidth', linWid)
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
ylim([-5, 5])
yline(1.96, '--', 'linewidth', linWid, 'color', 'red')
yline(-1.96, '--', 'linewidth', linWid, 'color', 'red')
xlim([-450, 3000])
export_fig([figSavePath 'Fig1ci_reactiveMethod_yes_mean' '.jpg'], '-r300')

figure('position', [0,0,1000,400])
hold on 
spread = 35;
tim = reactChan.HFB.encMulTim; 
for ii = 9:12
    hitIdx = find(reactChan.hits & reactChan.use); 
    RT = reactChan.encInfo(hitIdx(ii),4);
    plot(tim, allTrials(:,ii)+ii*spread, ...
        'color', purpleYellow((13-9)*60+1,:), 'linewidth', 2)
    scatter(repmat(RT,[50,1]), ii*linspace(spread-1.25,spread+1.25,50), 10, ...
        'red', 'filled')
    trial = allTrials(:,ii);            
    test = gausswin(11);
    test = test ./ sum(test); 
    trial = [zeros(5,1); trial; zeros(5,1)];
    smoothT = conv(trial, test, 'valid'); 
    coli = 14-ii;
    plot(tim, smoothT+ii*spread, ...
        'color', purpleYellow((1-1)*60+1,:), 'linewidth', 2)
    [maxVal, idx] = max(smoothT(find(tim==0):find(tim>=RT,1))); 
    testTim = tim(find(tim==0):find(tim>=RT,1));
    scatter(testTim(idx), smoothT(tim==testTim(idx))+ii*spread,...
            140,'k', 'filled' )
    scatter(testTim(idx), smoothT(tim==testTim(idx))+ii*spread,...
            100,[126, 165, 21]./255, 'filled' )

end
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca; ax.LineWidth=4;
xlim([-450, 3000])
yticks([])
export_fig([figSavePath 'Fig2a_HFBpeakDetection' '.jpg'], '-r300')


figure('position', [0,0,600,1200])
allTrials = [nonReact.HFB.subHit, nonReact.HFB.subMiss]; 
imagesc(nonReact.HFB.encMulTim, [], allTrials')
colormap(s2w2y)
clim([0,10])
colorbar
hold on 
xline(0, '--', 'linewidth', linWid, 'color', 'white')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
yticks([])
xlim([-450, 3000])
export_fig([figSavePath 'Fig1cii_reactiveMethod_non_heat' '.jpg'], '-r300')

figure('position', [0,0,600,400])
plot(nonReact.HFB.encMulTim, mean(allTrials,2), 'linewidth', linWid)
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
ylim([-5, 5])
yline(1.96, '--', 'linewidth', linWid, 'color', 'red')
yline(-1.96, '--', 'linewidth', linWid, 'color', 'red')
xlim([-450, 3000])
export_fig([figSavePath 'Fig1cii_reactiveMethod_non_mean' '.jpg'], '-r300')
    




