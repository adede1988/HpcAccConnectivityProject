function [] = makeLLplot(tim, hitVals, missVals, tVals,...
            clusterinfo, regSubs, path , roiLabs, memType)

 %set colors
colorNeg = {'#223C20', '#4C8D26', '#D5FB00', '#DE60CA', '#882380'};
colorPos = {'#100102', '#4B1E19', '#C0587E', '#FC8B5E', '#EA592A'};


f = figure('visible', false);
f.Position = [0 0 600 1300];
subplot 421
imagesc(tim, -150:150, squeeze(mean(hitVals)))
allmin = min(min(mean(hitVals), [], 'all'), min(mean(missVals),[], 'all')); 
allmax = max(max(mean(hitVals), [], 'all'), max(mean(missVals),[], 'all')); 
allmax = max([abs(allmin), allmax]); 
ylabel([roiLabs{1} ' leads      ' roiLabs{2} ' leads'])
title([memType ' Hit'])
yline(0)
xline(0)
caxis([-allmax, allmax])
colorbar

subplot 422
imagesc(tim, -150:150, squeeze(mean(missVals)))
ylabel([roiLabs{1} ' leads      ' roiLabs{2} ' leads'])
caxis([-allmax, allmax])
title([memType  ' Miss'])
yline(0)
xline(0)
colorbar

subplot 423
hold off
imagesc(tVals)
allmax = max([abs(min(tVals, [], 'all')), max(tVals, [], 'all') ]);
ylabel([roiLabs{1} ' leads      ' roiLabs{2} ' leads'])
caxis([-allmax, allmax])
title('t-value')
colorbar

cluSubi = 4;
Xidx = repmat(1:size(tVals,2), [size(tVals,1),1]);
Yidx = repmat(1:size(tVals,1), [size(tVals,2),1])';
negClu = []; 
if isfield(clusterinfo, 'neg_clusters')
    if isfield(clusterinfo.neg_clusters, 'p')
        negClu = find([clusterinfo.neg_clusters.p]<.05);
    end    
end
if ~isempty(negClu)
    for cc = 1:length(negClu)
        subplot 423
        tmp = clusterinfo.neg_clusters(negClu(cc)).inds;
        curp = ones(size(tmp)); 
        curp(tmp) = 0; 
        if cc < 6
            addRedOutline(curp', .05, colorNeg{cc});
        else
            addRedOutline(curp', .05, 'magenta');

        end


        if cluSubi < 9
            subplot(4,2,cluSubi)
            pairMeans = zeros((length(regSubs))*2,1); %hits, misses
            hmSort = pairMeans; 
            ti = 1; 
            for pi = 1:length(regSubs)
                pairVals = squeeze(hitVals(pi,:,:)); 
                pairMeans(ti) = mean(pairVals(tmp), 'all');
                hmSort(ti) = 1; 
                ti = ti+1; 
                pairVals = squeeze(missVals(pi,:,:)); 
                pairMeans(ti) = mean(pairVals(tmp), 'all');
                ti = ti+1; 
            end

            hold off

            b = boxchart(hmSort, pairMeans);
            b.MarkerStyle = 'none'; 
            xticks([0,1])
            xticklabels({'Miss', 'Hit'})
            meanX = round(mean(Xidx(tmp), 'all')); 
            meanY = round(mean(Yidx(tmp), 'all')); 
            LLvals = -150:150; 
            title(['time:  ' num2str(tim(meanX)) 'ms  leadLag: ' num2str(LLvals(meanY)) ' ms' ])
            PLH = pairMeans(hmSort==1); 
            PLM = pairMeans(hmSort==0); 
            randVals = (rand(length(PLH),1)-.5)*.5;
            hold on 
            scatter(randVals, PLM, 10,  'blue')
            scatter(randVals+1, PLH, 10, 'blue')
            
            for pi = 1:length(PLH)
                plot([0+randVals(pi),1+randVals(pi)], [PLM(pi),PLH(pi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
                
        
            end

            %make subject means
            subHits = zeros(length(unique(regSubs)),1); 
            subMisses = zeros(length(unique(regSubs)),1); 
            uniqueSubs = unique(regSubs); 
            for sub = 1:length(subHits)
                subHits(sub) = mean(PLH(regSubs==uniqueSubs(sub)));
                subMisses(sub) = mean(PLM(regSubs==uniqueSubs(sub)));
                plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

            end
            scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
            scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
            

            ylabel('leadLag strength (r)')



            ylabel('leadLag strength (r)')
            cluSubi = cluSubi + 1; 
        end
        
    end
end
posClu = []; 
if isfield(clusterinfo, 'pos_clusters')
    if isfield(clusterinfo.pos_clusters, 'p')
        posClu = find([clusterinfo.pos_clusters.p]<.05);
    end
end
if ~isempty(posClu)
    for cc = 1:length(posClu)
        subplot 423
        tmp = clusterinfo.pos_clusters(posClu(cc)).inds;
        curp = ones(size(tmp)); 
        curp(tmp) = 0; 
        if cc < 6
            addRedOutline(curp', .05, colorPos{cc});
        else
            addRedOutline(curp', .05, 'magenta');

        end

        if cluSubi < 9
            subplot(4,2,cluSubi)
            pairMeans = zeros((length(regSubs))*2,1); %hits, misses
            hmSort = pairMeans; 
            ti = 1; 
            for pi = 1:length(regSubs)
                pairVals = squeeze(hitVals(pi,:,:)); 
                pairMeans(ti) = mean(pairVals(tmp), 'all');
                hmSort(ti) = 1; 
                ti = ti+1; 
                pairVals = squeeze(missVals(pi,:,:)); 
                pairMeans(ti) = mean(pairVals(tmp), 'all');
                ti = ti+1; 
            end

            hold off

            b = boxchart(hmSort, pairMeans);
            b.MarkerStyle = 'none'; 
            xticks([0,1])
            xticklabels({'Miss', 'Hit'})
            meanX = round(mean(Xidx(tmp), 'all')); 
            meanY = round(mean(Yidx(tmp), 'all')); 
            LLvals = -150:150; 
            title(['time:  ' num2str(tim(meanX)) 'ms  leadLag: ' num2str(LLvals(meanY)) ' ms' ])
            PLH = pairMeans(hmSort==1); 
            PLM = pairMeans(hmSort==0); 
            randVals = (rand(length(PLH),1)-.5)*.5;
            hold on 
            scatter(randVals, PLM, 10,  'blue')
            scatter(randVals+1, PLH, 10, 'blue')
            
            for pi = 1:length(PLH)
                plot([0+randVals(pi),1+randVals(pi)], [PLM(pi),PLH(pi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
                
        
            end


            %make subject means
            subHits = zeros(length(unique(regSubs)),1); 
            subMisses = zeros(length(unique(regSubs)),1); 
            uniqueSubs = unique(regSubs); 
            for sub = 1:length(subHits)
                subHits(sub) = mean(PLH(regSubs==uniqueSubs(sub)));
                subMisses(sub) = mean(PLM(regSubs==uniqueSubs(sub)));
                plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

            end
            scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
            scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
            

            ylabel('leadLag strength (r)')
            cluSubi = cluSubi + 1; 
        end

    end
end
subplot 423
xticks([21,61,101])
xticklabels(tim([21,61,101]))
yticks([1:50:301])
yticklabels([-150:50:150])
yline(151)
xline(21)
 

export_fig(join([path memType '_' roiLabs{1} '_' roiLabs{2} '.jpg'],''), '-r300')
 
end