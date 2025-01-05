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




%% Figure 4A-D which frequencies have the most connectivity? and when? 
allConnections2 = load([datPre 'HFBConnections.mat']).allConnections2;
allConnections = load([datPre 'imageConnections.mat']).allConnections;
%connections relative to HFB by frequency: 
pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 1) ./ 41 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 1) ./ 41 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:4:20])
xticklabels(round(curDat.frex([1:4:20])))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
export_fig([figSavePath  'Fig4A_' 'CountSigConnections_byFrex_' ...
     '_HFBLocked.jpg'], '-r300')



%connections relative to image by frequency: 
pVals = squeeze(allConnections(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 1) ./ 139 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 1) ./ 139 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:4:20])
xticklabels(round(curDat.frex([1:4:20])))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;

export_fig([figSavePath 'Fig4C_' 'CountSigConnections_byFrex_' ...
     '_imageLocked.jpg'], '-r300')

%connections relative to HFB by time: 
pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
   
encSpect = sum(pVals(:,:,1), 2) ./ 20 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 2) ./ 20 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:10:41])
xticklabels([-500:250:500])
xline(21, 'linestyle', '--', 'linewidth', linWid, 'color', 'k')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;

export_fig([figSavePath 'Fig4B_' 'CountSigConnections_byTime' ...
     '_HFBLocked.jpg'], '-r300')

%connections relative to image by time: 
pVals = squeeze(allConnections(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 
x = curDat.enc_image_tim'; 
x(x<-450 | x>3000) = []; 
figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 2) ./ 20 ./ 25; 

y = encSpect; 
plot(x, y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 2) ./ 20 ./ 25;

y =retSpect;  
plot(x, y, 'color', 'blue', 'linewidth', 4)

xline(0, 'linestyle', '--', 'linewidth', linWid, 'color', 'k')
ylim([0, .08])
xlim([-450, 3000])
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;

export_fig([figSavePath 'Fig4D_' 'CountSigConnections_byTime' ...
     '_imageLocked.jpg'], '-r300')
















