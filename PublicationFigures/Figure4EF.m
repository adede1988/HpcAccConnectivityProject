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


%% Figure 4E and F schematic plots showing connections 

allConnections2 = load([datPre 'HFBConnections.mat']).allConnections2;
allConnections = load([datPre 'imageConnections.mat']).allConnections;

tim = curDat.enc_image_tim; 
tim(tim<-450 | tim>3000) = []; 
fn = [figSavePath...
        'Fig4F_schematic_image_']; 
fn2 = [figSavePath...
        'Fig4F_schematic_colKey.jpg']; 
fn3 = [figSavePath...
        'Fig4F_schematic_tKey.jpg']; 

%early encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-250 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc1.jpg'], ...
    fn2, fn3, datPre)

%mid encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=500 & tim<=1500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc2.jpg'], ...
    fn2, fn3, datPre)

%late encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=1500 & tim<=2500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc3.jpg'], ...
    fn2, fn3, datPre)

%early retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-250 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret1.jpg'], ...
    fn2, fn3, datPre)

%mid retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=500 & tim<=1500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret2.jpg'], ...
    fn2, fn3, datPre)

%late retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=1500 & tim<=2500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret3.jpg'], ...
    fn2, fn3, datPre)


%% schematics at HFB locked time

tim = [-500:25:500];  
fn = [figSavePath...
        'Fig4E_schematic_HFB_']; 
fn2 = [figSavePath...
        'Fig4E_schematic_colKey.jpg']; 
fn3 = [figSavePath...
        'Fig4E_schematic_tKey.jpg']; 
%early encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, true, [fn 'enc1.jpg'], ...
    fn2, fn3, datPre)

%mid encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc2.jpg'], ...
    fn2, fn3, datPre)

%late encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc3.jpg'], ...
    fn2, fn3, datPre)

%early retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret1.jpg'], ...
    fn2, fn3, datPre)

%mid retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret2.jpg'], ...
    fn2, fn3, datPre)

%late retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret3.jpg'], ...
    fn2, fn3, datPre)

%% same figures but with region colors
tim = [-500:25:500];  
fn = [figSavePath...
        'Fig4E_schematic_HFB_reg_']; 
fn2 = [figSavePath...
        'Fig4E_schematic_colKey.jpg']; 
fn3 = [figSavePath...
        'Fig4E_schematic_tKey.jpg']; 
%early encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc1.jpg'], ...
    fn2, fn3, datPre)

%mid encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc2.jpg'], ...
    fn2, fn3, datPre)

%late encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc3.jpg'], ...
    fn2, fn3, datPre)

%early retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret1.jpg'], ...
    fn2, fn3, datPre)

%mid retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret2.jpg'], ...
    fn2, fn3, datPre)

%late retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret3.jpg'], ...
    fn2, fn3, datPre)

















