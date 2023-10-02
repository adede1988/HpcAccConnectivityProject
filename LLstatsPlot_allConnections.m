

% LL stats plotting

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_LL_KEY_STATS';

statFiles = dir(path); 

statFiles(1:2) = []; 
LLdat = load([statFiles(1).folder '/' statFiles(1).name]).LLdat; 


%hit/miss, connection, offset, time
allConSub = zeros(2,20000,301, 141); 
% allConRet = zeros(2,20000,150, 121); 
%connection X summaryInfo
%1: subIDnum
%2: chanIDnum 1
%3: chanIDnum 2
%4: ROIi1
%5: ROIi2
allConSummary = zeros(20000, 5); 

ri = 1; 
for ii = 1:length(statFiles)
ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    if LLdat.reg1>LLdat.reg2
    hitVals = permute(squeeze(LLdat.regRes(1,:,:,:)), [3,1,2]); %only get the leading relationships
    missVals = permute(squeeze(LLdat.regRes(2,:,:,:)), [3,1,2]); %only get the leading relationships
    
    allConSub(1,ri:ri+size(hitVals,1)-1, :, :) = hitVals; 
    allConSub(2,ri:ri+size(missVals,1)-1,:,:) = missVals; 
    

    chani = cellfun(@(x) split(x, "_"), LLdat.chani, 'uniformoutput', false);
    chani1 = cellfun(@(x) str2num(x{1}), chani);
    chani2 = cellfun(@(x) str2num(x{2}), chani);

    allConSummary(ri:ri+LLdat.n_pair-1, 1) = LLdat.regSubs; 
    allConSummary(ri:ri+LLdat.n_pair-1, 2) = chani1; 
    allConSummary(ri:ri+LLdat.n_pair-1, 3) = chani2; 
    allConSummary(ri:ri+LLdat.n_pair-1, 4) = ones(LLdat.n_pair, 1) * LLdat.reg1; 
    allConSummary(ri:ri+LLdat.n_pair-1, 5) = ones(LLdat.n_pair, 1) * LLdat.reg2; 

    

%% retrieval! 

%     hitVals = permute(squeeze(LLdat.regRes2(1,:,:,:)), [3,1,2]);
%     missVals = permute(squeeze(LLdat.regRes2(2,:,:,:)), [3,1,2]);
%    
%     allConRet(1,ri:ri+size(hitVals,1)-1, :, :) = hitVals; 
%     allConRet(2,ri:ri+size(missVals,1)-1,:,:) = missVals; 

    ri = ri+size(hitVals,1); 
    end
end

%trim extra space in the output matrices 
allConSummary(ri:end,:) = []; 
allConSub(:,ri:end,:,:) = []; 

%cut off the trial ends, since these may be contaminated by flipping anyway
% allConSub(:,:,:,101:end) = []; 
tim = LLdat.encTim; %(1:100); 

%flatten into chanPair X offset/time
%eliminate within region pairs
distsHit = squeeze(allConSub(1, allConSummary(:,4)~=allConSummary(:,5), :, :)); 
distsHit = reshape(distsHit, [size(distsHit,1), prod(size(distsHit,[2,3])) ]); 
distsMiss = squeeze(allConSub(2, allConSummary(:,4)~=allConSummary(:,5), :, :)); 
distsMiss = reshape(distsMiss, [size(distsMiss,1), prod(size(distsMiss, [2,3]))]);
 
%get covariance
distsHitCov = distsHit*distsHit' ./ size(distsHit,1); 
distsMissCov = distsMiss*distsMiss' ./ size(distsMiss,1); 

% regularize R matrix
gamma = .01;
distsMissCov = distsMissCov *(1-gamma) + eye(length(distsMissCov))*gamma*mean(eig(distsMissCov));

%do GED
[V, L] = eig(distsHitCov, distsMissCov); 
[L, sidx] = sort(diag(L), 'descend'); 
V = V(:,sidx); 
tmpSum = allConSummary(allConSummary(:,4) ~= allConSummary(:,5), :);
cutVal = 95;

% https://academic.oup.com/bioinformatics/article/25/3/401/244239



%deal with sign uncertainty
for ii = 1:size(V,2)  
    w = V(:,ii); 
    [~, wi] = max(abs(w)); 
    if w(wi)<0
        w = w*-1; 
    end
    V(:,ii) = w; 
end

%find similarity of all connections to top 25 components (95% of variance)
cmp = (V(:,1:25)' * distsHit); %create the component "timeseries"
cmp = reshape(cmp, [size(cmp,1), size(allConSub, [3,4])] ); %reshape to offset X time
hitCons = squeeze(allConSub(1, allConSummary(:,4)~=allConSummary(:,5), :, :)); %get the hit connections
allSims = zeros(size(cmp,1), size(hitCons,1)); %initialize matrix of similarities (correlations)
parfor ii = 1:size(cmp,1)
    tic
    slice = allSims(ii,:); 
    curCMP = squeeze(cmp(ii,:,:)); 
    [coeff,score,latent,tsquared,explained,mu] = pca(curCMP);
    for jj = 1:size(hitCons,1)
        [coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(squeeze(hitCons(jj,:,:)));

        %% scratch work
%         figure
%         subplot 221
%         imagesc(curCMP)
%         subplot 222
%         imagesc(score(:,1) * coeff(:,1)')
%         subplot 223
%         imagesc(squeeze(hitCons(jj,:,:)))
%         subplot 224
%         imagesc(score2(:,1) * coeff2(:,1)')
        
%         matCor(score(:,1) * coeff(:,1)', score2(:,1) * coeff2(:,1)')


        %%
        PCAsims = zeros(10,1); 
        for cc = 1:10
            PCAsims(cc) = matCor(score(:,cc) * coeff(:,cc)', score2(:,cc) * coeff2(:,cc)');
        end
        slice(jj) = wmean(PCAsims, explained2(1:10)); 
%         slice(jj) = matCor(score(:,1) * coeff(:,1)', score2(:,1) * coeff2(:,1)');

%         matCor(curCMP, squeeze(hitCons(jj,:,:))); 
    end
    allSims(ii,:) = slice; 
    toc
end

%do empirical connections tend to have a singular component connection that
%they are similar to? 
%1: max cor val, 
%2: index of max cmp, 
%3: max - 2nd most, 
%4: max - mean, 
%5: 2nd most
maxSim = zeros(size(allSims,2), 5); 
for ii = 1:size(allSims,2)
    [maxSim(ii,1), maxSim(ii,2)] =  max(allSims(:,ii)); 
    [sortedCur, ~] = sort(allSims(:,ii), 'descend');
    maxSim(ii,3) = sortedCur(1) - sortedCur(2); 
    maxSim(ii,4) = allSims(maxSim(ii,2), ii) - mean(allSims(:,ii)); 
    maxSim(ii,5) = sortedCur(2); 

end

%% plotting 
for ii = 1:10
figure
% f = figure('visible', false);
% f.Position = [0 0 600 1000];
cmpidx = maxSim(:,1)>.3 & maxSim(:,2)==ii & maxSim(:,3)>.2;
cmp = mean(distsHit(cmpidx, :));
cmp = reshape(cmp, size(allConSub, [3,4]));
subplot(6, 2, [1,2])
imagesc(tim,-150:50:150, cmp)
propVar = L(ii) / sum(L); 
 ylabel(['reg2           reg1'])
title(['GED componenent: ' num2str(ii) '; sub HIT; prop Var:' num2str(round(propVar,2))])
yline(0)
allMax = max(cmp, [], 'all');
caxis([0, allMax])
colorbar



subplot(6,2,3)
hold off
histogram(maxSim(cmpidx, 1), [0:.05:1])
hold on 
histogram(maxSim(cmpidx, 5), [0:.05:1])

subplot(6,2,4)
curSum = tmpSum(cmpidx,:); 
propCons = zeros(11); 
for r1 = 1:11
    for r2 = 1:11
        nom = sum(curSum(:,4) == r1 & curSum(:,5) == r2);
        den = sum(tmpSum(:,4) == r1 & tmpSum(:,5) == r2);
        propCons(r1, r2) = nom/den; 

    end
end

% propCons(3,1) = 10; 
imagesc(propCons)
ROInames = {LLdat.aggTargs.ROI};
xticks(1:11)
yticks(1:11)
xticklabels(ROInames)
yticklabels(ROInames)
ylabel("reg1")
xlabel("reg2")
caxis([sum(cmpidx)/length(cmpidx), max(propCons, [], 'all')])
colorbar

%get the mean of the top Miss empirical maps
cmp = mean(distsMiss(cmpidx, :));
cmp = reshape(cmp, size(allConSub, [3,4])); 
subplot(6,2,[5,6])
imagesc(tim,-150:50:150, cmp)
caxis([0, allMax])
colorbar
ylabel(['reg2           reg1'])
title(['sub MISS'])


%take three random samples to fill in the bottom three rows of the figure
testi = randsample(find(cmpidx), 3);

subplot(6,2,[7,8])
curCmp = reshape(distsHit(testi(1), :),size(allConSub, [3,4])) ;
imagesc(tim,-150:50:150, curCmp)

subplot(6,2,[9,10])
curCmp = reshape(distsHit(testi(2), :),size(allConSub, [3,4])) ;
imagesc(tim,-150:50:150, curCmp)

subplot(6,2,[11,12])
curCmp = reshape(distsHit(testi(3), :),size(allConSub, [3,4])) ;
imagesc(tim,-150:50:150, curCmp)




end



test = squeeze(mean(distsHit(hit>prctile(hit, cutVal), :),1)); 
test2 = squeeze(mean(distsMiss(hit>prctile(hit, cutVal), :),1)); 
test = reshape(test, size(allConSub, [3,4])); 
test2 = reshape(test2, size(allConSub, [3,4])); 
subplot(4, 2, [1,2])

imagesc(tim, -150:50:150, test)
propVar = L(ii) / sum(L); 
 ylabel(['reg2 leads           reg1 leads'])
title(['GED componenent: ' num2str(ii) '; sub HIT; prop Var:' num2str(round(propVar,2))])
allMax = max([max(test,[], 'all'), max(test2,[], 'all')]);
caxis([.02, allMax*.9])
colorbar

subplot(4,2,[3,4])
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, test2)
propVar = L(ii) / sum(L); 
 ylabel(['reg2 leads           reg1 leads'])
title(['GED componenent: ' num2str(ii) '; sub MISS; prop Var:' num2str(round(propVar,2))])
caxis([.02, allMax*.9])
colorbar


subplot(4,2,[5,6])
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, test - test2)
propVar = L(ii) / sum(L); 
 ylabel(['reg2 leads           reg1 leads'])
title(['difference'])
colorbar


subplot(4,2,7)

cmpidx = hit>prctile(hit,cutVal); 
curSum = tmpSum(cmpidx,:); 
propCons = zeros(11); 
for r1 = 1:11
    for r2 = 1:11
        nom = sum(curSum(:,4) == r1 & curSum(:,5) == r2);
        den = sum(tmpSum(:,4) == r1 & tmpSum(:,5) == r2);
        propCons(r1, r2) = nom/den; 

    end
end

% propCons(3,1) = 10; 
imagesc(propCons)
ROInames = {LLdat.aggTargs.ROI};
xticks(1:11)
yticks(1:11)
xticklabels(ROInames)
yticklabels(ROInames)
ylabel("reg1")
xlabel("reg2")
colorbar

% subplot(4,2,8)
% hitSet = squeeze(distsHit(hit>prctile(hit, 90), :));
% missSet = squeeze(distsMiss(hit>prctile(hit, 90), :)); 
% tVals = myArrayT(hitSet, missSet, 1); 
% test = max(tVals);


export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\GEDFigs\cmp' num2str(ii) '.jpg'] )


end


for ii = 1:10
f = figure('visible', false);
f.Position = [0 0 600 1000];
hit = V(:,ii)' * distsHitCov; 
test = squeeze(mean(distsHit(hit>prctile(hit, cutVal), :),1)); 
test2 = squeeze(mean(distsMiss(hit>prctile(hit, cutVal), :),1)); 
test = reshape(test, size(allConSub, [3,4])); 
test2 = reshape(test2, size(allConSub, [3,4])); 
testSum = tmpSum(hit>prctile(hit,cutVal),:);
subplot(5, 1, 1)
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, test)
propVar = L(ii) / sum(L); 
 ylabel(['reg2 leads           reg1 leads'])
title(['GED componenent: ' num2str(ii) '; sub HIT; prop Var:' num2str(round(propVar,2))])
allMax = max([max(test,[], 'all'), max(test2,[], 'all')]);
caxis([.02, allMax*.9])
colorbar


topidx = find(hit>prctile(hit, cutVal)); 
sampidx = randsample(1:size(testSum,1), 10); 

for si = 1:10

test = distsHit(topidx(sampidx(si)), :); 
test = reshape(test, size(allConSub, [3,4])); 
subplot(5,2,si)
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, test)
propVar = L(ii) / sum(L); 
 ylabel([ROInames{testSum(sampidx(si), 4)} '            ' ROInames{testSum(sampidx(si), 5)} ])
 title(testSum(sampidx(si), 1))


end
% 
% test = distsHit(topidx(sampidx(2)), :); 
% test = reshape(test, size(allConSub, [3,4])); 
% subplot(5,2,[3,4])
% % hit = (V(:,ii)' * distsHit); 
% % hit = reshape(hit, size(allConSub, [3,4])); 
% imagesc(LLdat.encTim, -150:50:150, test)
% propVar = L(ii) / sum(L); 
%  ylabel([ROInames{testSum(sampidx(1), 4)} '            ' ROInames{testSum(sampidx(2), 5)} ])
% colorbar


export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\GEDFigs\cmp' num2str(ii) '_examples' '.jpg'] )


end

%use generalized eigen decomposi
%finding the most similar 
for ii = 1:10
f = figure('visible', false);
f.Position = [0 0 600 1000];
hit = V(:,ii)' * distsHitCov; 
test = squeeze(mean(distsHit(hit>prctile(hit, cutVal), :),1)); 
test2 = squeeze(mean(distsMiss(hit>prctile(hit, cutVal), :),1)); 
test = reshape(test, size(allConSub, [3,4])); 
test2 = reshape(test2, size(allConSub, [3,4])); 
subplot(5, 2, [1:2])
%get scaled component values
cmp = (V(:,ii)' * distsHit); 
cmp = cmp - min(cmp, [], 'all');
cmp = cmp ./ max(cmp, [], 'all'); 
%make scaled version of all Hit maps
distsHitScaled = bsxfun(@minus, distsHit,  min(distsHit,[], 2)); 
distsHitScaled = distsHitScaled ./ max(distsHitScaled, [], 2); 

diffs = sqrt(sum((bsxfun(@minus, distsHit, cmp)).^2,2)); 
cmp = reshape(cmp, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, cmp)
propVar = L(ii) / sum(L); 
 ylabel(['reg2 leads           reg1 leads'])
title(['GED componenent: ' num2str(ii) '; sub HIT; prop Var:' num2str(round(propVar,2))])
% allMax = max([max(test,[], 'all'), max(test2,[], 'all')]);
% caxis([.02, allMax])
subplot(5,2,3)
histogram(diffs)

[difSorted, order] = sort(diffs); %what are the indices of the most similar maps
top10idx = order(1:round(length(order)/40));
top10idx2 = find(hit>prctile(hit,97.5));
xline(difSorted(round(length(order)/40)))

subplot(5,2,4)
idxAgree = arrayfun(@(x) ismember(x, top10idx), top10idx2); %agreement is low! 
%use the distance metric rather than the top weights! 

curSum = tmpSum(top10idx,:); 
propCons = zeros(11); 
for r1 = 1:11
    for r2 = 1:11
        nom = sum(curSum(:,4) == r1 & curSum(:,5) == r2);
        den = sum(tmpSum(:,4) == r1 & tmpSum(:,5) == r2);
        propCons(r1, r2) = nom/den; 

    end
end

% propCons(3,1) = 10; 
imagesc(propCons)
ROInames = {LLdat.aggTargs.ROI};
xticks(1:11)
yticks(1:11)
xticklabels(ROInames)
yticklabels(ROInames)
ylabel("reg1")
xlabel("reg2")
colorbar


testSum = tmpSum(top10idx,:);
topidx = top10idx; 
sampidx = randsample(1:size(testSum,1), 10); 

for si = 1:3

test = distsHit(topidx(si), :); 
test = reshape(test, size(allConSub, [3,4])); 
subplot(5,2,[5+(si-1)*2, 5+(si-1)*2+1])
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(tim, -150:50:150, test)
propVar = L(ii) / sum(L); 
 ylabel([ROInames{testSum(si, 4)} '            ' ROInames{testSum(si, 5)} ])
 title(testSum(si, 1))


end


export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\GEDFigs\cmp' num2str(ii) '.jpg'] )


end




%use generalized eigen decomposition to ask where the leadLag energy is that differentiates between hits and misses

% 
% [coeff,score,latent,tsquared,explained,mu] = pca(distsHit);
% 
% 
% % distsHit = corr(distsHit); 
% test = DBscanDynamicEpi(distsHit, 10, 3, 1, 1); 

%how do pairs of electrodes covary with each other across the offset X time
%space





