

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'mni2atlas']))
% addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
% ft_defaults
datFolder = [datPre 'HFB_singleTrial\out']; 
regFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({regFiles.name}, '.mat'));
regFiles = regFiles(test); 


for ii = 1: length(regFiles)
    fileBits = split(regFiles(ii).name, 'stat'); 
    regFiles(ii).stat = str2num(fileBits{2}(1));
    fileBits = split(fileBits{2}, '_'); 
    regFiles(ii).reg = fileBits{2}; 
    fileBits = split(fileBits{3}, '.'); 
    regFiles(ii).phase = fileBits{1}; 


end


regions = unique({regFiles.reg});
phases = {'sub', 'ret'}; 
for reg = 1:length(regions)

    for ph = 1:2
        try
        curFiles = regFiles(ismember({regFiles.reg}, regions{reg}) & ...
                            ismember({regFiles.phase}, phases{ph})); 
        HFBdat = load([curFiles(1).folder '/' curFiles(1).name]).statInfo; 
        tmp = load([curFiles(2).folder '/' curFiles(2).name]).statInfo;

        if isfield(tmp, 'tVals_image')
            HFBdat.tVals_image = tmp.tVals_image; 
            HFBdat.hitVals_image = tmp.hitVals_image; 
            HFBdat.missVals_image = tmp.missVals_image; 
            HFBdat.p_image = tmp.p_image; 
        else
            HFBdat.tVals_HFB = tmp.tVals_HFB; 
            HFBdat.hitVals_HFB = tmp.hitVals_HFB; 
            HFBdat.missVals_HFB = tmp.missVals_HFB; 
            HFBdat.p_HFB = tmp.p_HFB; 


        end

        tim = HFBdat.tim; 
         hits = HFBdat.hitVals_image; 
        misses = HFBdat.missVals_image;
      
      

        tim(tim<-450 | tim>3000) = [];
        latency = HFBdat.hitLat; 
        latency2 = HFBdat.missLat; 
    figure('visible', true, 'position', [0,0,600,1000])

    subplot(4,2,5)
    x = tim(10:end-10); 
    y = hits(10:end-10);
%     e = std(hits(10:end-10,:),[],2) ./ sqrt(length(latency));
    plot(x, y, 'color', 'blue', 'linewidth', 2); 
    hold on 

    y = misses(10:end-10);
    plot(x,y, 'color', 'red', 'linewidth', 2); 
    
    xlim([tim(10), tim(end-10)])

    subplot(4,2,7)

    y1 = hits(10:end-10);
    plot(x, y1 - y, 'color', 'k', 'linewidth', 2); 
    ylim([-max(abs(y1-y)), max(abs(y1-y))])
    yline(0, '--')
    hold on 
    x(HFBdat.p_image(10:end-10)>.05) = []; 

    plot(x, zeros(length(x),1), 'color', 'red', 'linewidth', 5)

    xlim([tim(10), tim(end-10)])


    subplot(4, 2, [1,3])
    [~, order] = sort(latency); 
    hits = HFBdat.hits;
    misses = HFBdat.misses;
    hits(tim<-450 | tim>3000,:) = []; 
    misses(tim<-450 | tim>3000,:) = []; 

    imagesc(tim(10:end-10), [], hits(10:end-10,order)')
    caxis([-10,10])
    hold on 
    xline(0, 'color', 'green', 'linewidth', 3)

    title([regions{reg} ' ' phases{ph} ' Hit'])
  

    subplot(4, 2, [2,4])


    [~, order] = sort(latency2);
    imagesc(tim(10:end-10), [], misses(10:end-10,order)')
       caxis([-10,10])
    xline(0, 'color', 'green', 'linewidth', 3)
    title([regions{reg} ' ' phases{ph} ' Miss'])

    subplot(4,2,6)
     hold off
     w = 20; 

     y = HFBdat.hitVals_HFB; 
     x = [-w*25:25:w*25]';
     plot(x, y, 'color', 'blue', 'linewidth', 2); 
     hold on 

     y2 = HFBdat.missVals_HFB; 
     x = [-w*25:25:w*25]';
     plot(x, y2, 'color', 'red', 'linewidth', 2); 

     subplot(4,2,8)
      
    plot(x, y - y2, 'color', 'k', 'linewidth', 2); 
    ylim([-max(abs(y - y2)), max(abs(y - y2))])
    yline(0, '--')
    hold on 
    x(HFBdat.p_HFB>.05) = []; 

    plot(x, zeros(length(x),1), 'color', 'red', 'linewidth', 5)

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\singleTrialStats\' ...
        regions{reg} '_' phases{ph} '.jpg'], '-r300')
    close all
        catch
        end


    end
end
















