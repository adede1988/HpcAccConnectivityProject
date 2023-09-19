

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

allSigEnc = zeros(length(LLdat.encTim), 11); 
allSigRet = zeros(length(LLdat.retTim), 11); 
regNames = cell(11,1); 

%hit/miss X region X early/late X encode/retrieve
%early = 0-1000ms; late = 1000-2000ms
allMeanEarlyLate = zeros(2, 11, 2, 2); 

allConSub = zeros(2,10000,301, 141); 
allConRet = zeros(2,10000,301, 121); 
ri = 1; 
for ii = 1:length(statFiles)
ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    hitVals = permute(squeeze(LLdat.regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(LLdat.regRes(2,:,:,:)), [3,1,2]);
    
    allConSub(1,ri:ri+size(hitVals,1)-1, :, :) = hitVals; 
    allConSub(2,ri:ri+size(missVals,1)-1,:,:) = missVals; 
  


%% retrieval! 

    hitVals = permute(squeeze(LLdat.regRes2(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(LLdat.regRes2(2,:,:,:)), [3,1,2]);
   
    allConRet(1,ri:ri+size(hitVals,1)-1, :, :) = hitVals; 
    allConRet(2,ri:ri+size(missVals,1)-1,:,:) = missVals; 

    ri = ri+size(hitVals,1); 
end


distsHit = squeeze(allConSub(1, :, :, :)); 
distsHit = reshape(distsHit, [size(distsHit,1), prod(size(distsHit,[2,3])) ]); 
distsMiss = squeeze(allConSub(2, :, :, :)); 
distsMiss = reshape(distsMiss, [size(distsMiss,1), prod(size(distsMiss, [2,3]))]);
 
distsHitCov = cov(distsHit'); 
distsMissCov = cov(distsMiss'); 
[V, D] = eig(distsHitCov, distsMissCov); 
[L, sidx] = sort(diag(D), 'descend'); 
V = V(:,sidx); 

for ii = 1:10
figure
subplot 211
hit = (V(:,ii)' * distsHit); 
hit = reshape(hit, size(allConSub, [3,4])); 
imagesc(LLdat.encTim, [], hit)
propVar = L(ii) / sum(L); 
title(['GED componenent: ' num2str(ii) ' prop Var:' num2str(round(propVar,2))])
colorbar

hit = V(:,ii)' * distsHitCov; 
test = squeeze(mean(distsHit(hit>prctile(hit, 90), :),1)); 
test = reshape(test, size(allConSub, [3,4])); 
subplot 212
imagesc(test)
colorbar

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





