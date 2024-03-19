function [] = connectionGraph(curDat, ...
                regions, phase, stat, cc, pp, datPre)

reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions));

confrex = logspace(log10(2), log10(25), 20);
stds = logspace(log10(3),log10(5),20)./(2*pi*confrex);
chanFiles = dir([datPre 'CHANDAT/finished/']);



%what's the mean frequency of this connection? 
weight = sum(curDat.([phase '_' stat '_clust']).pos_clusters(cc).inds);
meanFrex = sum((confrex .* weight)) / sum(weight);
[~, freqi] = min(abs(confrex - meanFrex));


%what's the mean time index for this connection? 
weight = sum(curDat.([phase '_' stat '_clust']).pos_clusters(cc).inds'); 
ti = round(sum([1:length(weight)] .* weight) / sum(weight)); 




    
    



%loop pairs to calculate AUC for each channel pair
if curDat.chiVals(pp,1) < 10
    chi = ['00' num2str(curDat.chiVals(pp,1))]; 
elseif curDat.chiVals(pp,1) < 100
    chi = ['0' num2str(curDat.chiVals(pp,1))]; 
else
    chi = num2str(curDat.chiVals(pp,1)); 
end
ch1 = load([datPre 'CHANDAT/finished/chanDat_' ...
    curDat.subVals{pp} '_' chi '.mat']).chanDat; 


%some cases are flipped because in
%connectionIntegrate.m line 25 I had to flip regions
%sometimes for the image locked data since timing
%doesn't matter for these, so I avoided repeat
%calculation by flipping for half of the pairs, but I
%then wrote down the flipped chi values for these and
%it does matter here where I'm reconstructing file
%names for HFB locked connectivity
if ~strcmp(ch1.labels{curDat.chiVals(pp,1),3}, regions{reg1})
    if curDat.chiVals(pp,2) < 10
        chi = ['00' num2str(curDat.chiVals(pp,2))]; 
    elseif curDat.chiVals(pp,2) < 100
        chi = ['0' num2str(curDat.chiVals(pp,2))]; 
    else
        chi = num2str(curDat.chiVals(pp,2)); 
    end
    ch1 = load([datPre 'CHANDAT/finished/chanDat_' ...
        curDat.subVals{pp} '_' chi '.mat']).chanDat;
    tmp = curDat.chiVals(pp,1); 
    curDat.chiVals(pp,1) = curDat.chiVals(pp,2); 
    curDat.chiVals(pp,2) = tmp; 
end


idx = cellfun(@(x) ~isempty(strfind(x, ch1.subID)), ...
        {chanFiles.name}); 
subFiles = chanFiles(idx); 
if strcmp(stat, 'HFB')
    HFBflag = true;
else
    HFBflag = false; 
end
%chan X chan matrix of ppc connection metrics at time and frequency
%of target connection
%undirected connectivity
%for HFB connections, connections will have to be calcualted from
%scratch since if current channel is X, but I want connections
%between Y and Z at the time of X's HFB peak, then that won't be
%possible without a fresh calculation from scratch. 
hitMat = zeros(length(subFiles)); 
missMat = zeros(length(subFiles));
disp(['channel count: ' num2str(length(subFiles))])
for ii = 1:length(subFiles)
    
    
    chii = load([subFiles(ii).folder '/' ...
                 subFiles(ii).name]).chanDat;
    for jj = 1:length(subFiles)
        if ii>jj %don't repeat work! 
          
        

    if HFBflag
    %PPC values need to be calcualted for HFB aligned data
        chjj = load([subFiles(jj).folder '/' ...
                         subFiles(jj).name]).chanDat;
    
        if strcmp(phase, 'enc')
            dat1 = chii.enc; 
            dat2 = chjj.enc; 
            di = chii.ISPC.encdi; 
            hitTrials = ch1.use & ch1.hits;
            missTrials = ch1.use & ch1.misses;
            HFBLatHit = ch1.HFB_lat.subHit;
            HFBLatMiss = ch1.HFB_lat.subMiss;
            outTim = ch1.HFB.encMulTim;
             
        else
            dat1 = chii.retOn; 
            dat2 = chjj.retOn; 
            di = chii.ISPC.ondi; 
            hitTrials = ch1.retInfo(:,1)==1;
            missTrials = ch1.retInfo(:,1)==2;
            HFBLatHit = ch1.HFB_lat.retHit;
            HFBLatMiss = ch1.HFB_lat.retMiss;
            outTim = ch1.HFB.onMulTim;
    
        end
        [~, ppcHit] = getClustAngDat(dat1, dat2, confrex, ...
            20, stds,...
            1000, di, hitTrials, ...
            HFBLatHit, ...
            outTim, HFBflag, freqi, ti);
        hitMat(ii,jj) = ppcHit; 
        
        [~, ppcMiss] = getClustAngDat(dat1, dat2, confrex, ...
            20, stds,...
            1000, di, missTrials, ...
            HFBLatMiss, ...
            outTim, HFBflag, freqi, ti);
        missMat(ii,jj) = ppcMiss; 
    else
        outTim = ch1.HFB.encMulTim;
        %PPC values can be pulled from prior calculations since this is
        %imaged locked time
        if strcmp(phase, 'enc')
           
            hitMat(ii,jj) = chii.ISPC.subHit(jj, ti, freqi, 2);
            missMat(ii,jj) = chii.ISPC.subMiss(jj, ti, freqi, 2);

        else

            hitMat(ii,jj) = chii.ISPC.hit_on(jj, ti, freqi, 2);
            missMat(ii,jj) = chii.ISPC.miss_on(jj, ti, freqi, 2);

        end
            

    end


        end
    end
    disp([num2str(ii) ': ' num2str(round(toc)/60)])
end
hitMat = hitMat + hitMat'; 
missMat = missMat + missMat'; 
N = size(hitMat, 1); 

test = hitMat; 
test(test <0) = 0; 
test = sqrt(test); 
test = 1./test; 

BC_hit=betweenness_wei(test) ./ ((N-1)*(N-2));
ST_hit=strengths_und(hitMat) ./ N; 

test = missMat; 
test(test <0) = 0; 
test = sqrt(test); 
test = 1./test; 

BC_miss=betweenness_wei(test) ./ ((N-1)*(N-2));
ST_miss=strengths_und(missMat) ./ size(missMat,1); 



outDat = struct; 

if HFBflag
    tim = [-500:25:500]; 
else
    tim = outTim; 
end
outDat.subID = curDat.subVals{pp}; %participant ID
outDat.chi = curDat.chiVals(pp,1); %originating channel index
outDat.chi2 = curDat.chiVals(pp,2); %partner channel index
outDat.HFB_Image = stat; %HFB / Image timing
outDat.encRet = phase; %enc / ret
outDat.reg1 = regions{reg1}; 
outDat.reg2 = regions{reg2}; 
outDat.labels = ch1.labels(:,3); %will be useful for sorting later
outDat.hitBC = BC_hit; %BC for all channels on hit trials
outDat.missBC = BC_miss; %BC for all channels on miss trials
outDat.hitST = ST_hit; %ST for all channels on hit trials
outDat.missST = ST_miss; %ST for all channels on miss trials
outDat.tVal = curDat.([phase '_' stat '_tVal'])(ti,freqi);
outDat.time = tim(ti); % time in ms 
outDat.freq = meanFrex; % mean frequency (HFB uses 120)
outDat.hitMat = hitMat; 
outDat.missMat = missMat; 

save([datPre 'graphAnalysis/out/' regions{reg1} '_' regions{reg2} '_'...
    phase '_' stat '_' num2str(cc) '_' num2str(pp) '.mat'], 'outDat')



end