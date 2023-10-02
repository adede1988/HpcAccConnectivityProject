
%TF stats plot
clear
codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\TF_KEY_STATS';

statFiles = dir(path); 

statFiles(1:2) = []; 
frex = logspace(log10(2),log10(80),100);


for ii = 1:length(statFiles)
    ii
    TFdat = load([statFiles(ii).folder '/' statFiles(ii).name]).TFdat; 
    roi_idx = TFdat.reg1; 
    roiLabs = TFdat.aggTargs(roi_idx).ROI;     
  
    f = figure('visible', false);
    f.Position = [0 0 1000 600];

    hitVals = permute(squeeze(TFdat.regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(TFdat.regRes(2,:,:,:)), [3,1,2]);
    tim = TFdat.encTim; 

    allmax = max([max(mean(hitVals,1), [], 'all'), max(mean(missVals,1), [], 'all')]);
    allmin = min([min(mean(hitVals,1), [], 'all'), min(mean(missVals,1), [], 'all')]);
    allmax = max([allmax, abs(allmin)]);
    subplot 321
    hold off
    imagesc(tim, [], squeeze(mean(hitVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([roiLabs ': subsequent hit z-score power'])
    colorbar
    caxis([-allmax, allmax])


    subplot 323
    hold off
    imagesc(tim, [], squeeze(mean(missVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([roiLabs ': subsequent miss z-score power'])
    colorbar
    caxis([-allmax, allmax])

    subplot 325
    hold off
    imagesc(TFdat.tVals_sub')
    test = addRedOutline(TFdat.p_sub, .2, 'green');
    test = addRedOutline(TFdat.p_sub, .05, 'red');
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    xticks([1:20:141])
    xticklabels(tim([1:20:141]))
    title([roiLabs ': sub hit v. miss t-values'])
    colorbar
  

    hitVals = permute(squeeze(TFdat.regRes2(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(TFdat.regRes2(2,:,:,:)), [3,1,2]);
    tim = TFdat.retTim; 

    allmax = max([max(mean(hitVals,1), [], 'all'), max(mean(missVals,1), [], 'all')]);
    allmin = min([min(mean(hitVals,1), [], 'all'), min(mean(missVals,1), [], 'all')]);
    allmax = max([allmax, abs(allmin)]);
    subplot 322
    hold off
    imagesc(tim, [], squeeze(mean(hitVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([roiLabs ': retrieval hit z-score power'])
    colorbar
    caxis([-allmax, allmax])


    subplot 324
    hold off
    imagesc(tim, [], squeeze(mean(missVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([roiLabs ': retrieval miss z-score power'])
    colorbar
    caxis([-allmax, allmax])

    subplot 326
    hold off
    imagesc(TFdat.tVals_ret')
    test = addRedOutline(TFdat.p_ret, .2, 'green');
    test = addRedOutline(TFdat.p_ret, .05, 'red');
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    xticks([1:20:121])
    xticklabels(tim([1:20:121]))
    title([roiLabs ': ret hit v. miss t-values'])
    colorbar




    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' roiLabs '.jpg'], '-r300')

end





