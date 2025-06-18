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
addpath([codePre 'fieldtrip-20230118'])
ft_defaults;
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


%% Figure 5A HIP ACC example connection



hip_enc_HFB_file = 'hip_acc_ret_HFB_2_21.mat';


exampDat = load([datPre...
                 hip_enc_HFB_file]).outDat;

missName = ['Fig5A_miss_'  'hip_acc_ret_HFB'];
hitName = ['Fig5A_hit_'  'hip_acc_ret_HFB']; 

hipModel =  stlread('G:\My Drive\GitHub\HpcAccConnectivityProject\bilatHippo.stl');


set(0, 'defaultfigurewindowstyle', 'normal')


ftpath = [codePre 'fieldtrip-20230118']; 
yeo7 = ft_read_atlas(...
    fullfile(ftpath, ...
    'template/atlas/yeo/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'));
hold on 
test = yeo7.tissue>0; 
indices = zeros(1003785, 3); 
vals = zeros(1003785, 1); 
indi = 1; 
for nn = 1:256
    for jj = 1:256
        for kk = 1:256
            if test(nn,jj,kk)
                indices(indi, :) = [nn*-1+127,...
                                    kk*-1+145,...
                                    jj*-1+147]; 
                vals(indi) = yeo7.tissue(nn, jj, kk); 
                indi = indi + 1; 
            end
        end
    end
end
figure('visible', true, 'position', [50,50,1000,1000])
test = randsample(1:indi, 10000, false); 
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)






view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

hitBin = exampDat.hitMat>.10; 
missBin = exampDat.missMat>.10; 
disp(sum(hitBin, 'all'))
sum(missBin, 'all')

for nn = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if hitBin(nn,jj)
            plot3([exampDat.chanpos(nn,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(nn,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(nn,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end

export_fig([figSavePath  hitName '.jpg'], '-r300')


figure('visible', true, 'position', [50,50,1000,1000])
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)

for nn = 1:size(missBin,1)
    for jj = 1:size(missBin,1)
        if missBin(nn,jj)
            plot3([exampDat.chanpos(nn,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(nn,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(nn,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])

set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

export_fig([figSavePath  missName '.jpg'], '-r300')


%% same data but image locked and averaged over 1.5 seconds

missName = ['Fig5A_miss_'  'hip_acc_ret_IMAGE'];
hitName = ['Fig5A_hit_'  'hip_acc_ret_IMAGE']; 

imHit = exampDat.imHitMat; 
imMiss = exampDat.imMissMat; 


figure('position', [50,50,1000,1000])
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')


view([250,-50, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 


hitBin = imHit>.0686; 
missBin = imMiss>.0686; 
sum(hitBin, 'all')
sum(missBin, 'all')
for ii = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if hitBin(ii,jj)
            plot3([exampDat.chanpos(ii,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(ii,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(ii,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
export_fig([figSavePath  hitName '.jpg'], '-r300')

figure('position', [50,50,1000,1000])
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')


view([250,-50, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 


for ii = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if missBin(ii,jj)
            plot3([exampDat.chanpos(ii,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(ii,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(ii,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
export_fig([figSavePath  missName '.jpg'], '-r300')

%% Supplemental Figure 10B histogram of connection strengths 

hitVal = triu(imHit,1); 
missVal = triu(imMiss,1);
hitVal = hitVal(:); 
missVal = missVal(:); 
hitVal(hitVal==0) = []; 
missVal(missVal==0) = []; 


figure
histogram(hitVal, [-.1:.01:.2], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.01:.2], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;

export_fig([figSavePath  'supFig10B_mainConnectionHist' '.jpg'], '-r300')

figure
histogram(hitVal, [-.1:.02:.7], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.02:.7], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
xlim([.025, .7])
ylim([0,10])
export_fig([figSavePath  'supFig10B_detailZoomConnectionHist' '.jpg'], '-r300')

figure
histogram(hitVal, [-.1:.01:.7], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.01:.7], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
xlim([-.025, .025])













