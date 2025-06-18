% channel HFB assessment


%github access token: 
%  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
%  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


%local paths: 


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'mni2atlas']))
% addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
% ft_defaults
datFolder = [datPre 'CHANDAT/finished']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 







%% loop over channels and get single trial level data for HFB and TF

regions = {'acc', 'dlPFC', 'hip', ...
    'lTemp', 'iTemp', 'mtl', 'pcc', 'pPFC', 'vis'}; 

regCounts = zeros(length(chanFiles), 1); 
%regions are encoded in different powers of 10
%10^regidx + code
%codes for reactivity: 
%1: reactive at encoding
%2: reactive at retrieval
%3: reactive at both
%0: reactive at neither

parfor ii = 1:length(chanFiles)
ii
chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
lab = chanDat.labels{chanDat.chi,3}; 

T = sum(chanDat.retInfo(:,1)==1 | chanDat.retInfo(:,1)==2); 
Hr = sum(chanDat.retInfo(:,1)==1) / T; 
T = sum(chanDat.retInfo(:,1)==3 | chanDat.retInfo(:,1)==4); 
F = sum(chanDat.retInfo(:,1)==4);
if F == 0
    acc =  Hr - F/T; 
    F = 1; 
else
    acc =  Hr - F/T; 
end
Fr = F / T; 

d = norminv(Hr) - norminv(Fr); 

if acc>0 && chanDat.age > 13 && ~strcmp(lab, 'ZZZ')%memory and age 
    regidx = find(cellfun(@(x) strcmp(x, lab), regions));
    if sum(chanDat.reactiveRes(1:2))>0 %reactive during encoding
        regCounts(ii) = 10.^regidx + 1; 
    end
    if sum(chanDat.reactiveRes(3:4))>0 %reactive during retrieval 
        regCounts(ii) = 10.^regidx + 2; 
    end

    if sum(chanDat.reactiveRes(1:2))>0 && ... 
        sum(chanDat.reactiveRes(3:4))>0 %reactive both 
        regCounts(ii) = 10.^regidx + 3; 
    end
    if sum(chanDat.reactiveRes(1:4))==0
        regCounts(ii) = 10.^regidx; 
    end
end
end

regCountsExtract = zeros(9, length(regCounts)); 
for ii = 1:length(regCounts)

    if regCounts(ii) > 0
    code = mod(regCounts(ii), 10);
    reg = regCounts(ii) - code; 
    reg = log10(reg);
    regCountsExtract(reg, ii) = code+1; 
    %convert code: 
    %1: non react
    %2: encode react
    %3: retrieve react
    %4: both react

    else
        regCountsExtract(:,ii) = -1; 
    end
end

%reg, encRet, type, count
regCountsExtract(:, sum(regCountsExtract)<0) = []; 
n = 36;
aovDat = table;
aovDat.encRet = repmat("askj", n,1); 
aovDat.reg = repmat("askj", n,1); 
aovDat.count = zeros(n,1);
ai = 1; 

for reg = 1:9
    %non = 1
    %enc = 2
    %ret = 3
    %both = 4
    curReg = regCountsExtract(reg,:);
    curReg(curReg==0) = [];
    aovDat.encRet(ai) = 'encReact';  
    aovDat.count(ai)= sum(curReg==2);
    aovDat.reg(ai) = regions{reg}; 
    ai = ai+1; 

    aovDat.encRet(ai) = 'retReact'; 
    aovDat.count(ai)= sum(curReg==3);
    aovDat.reg(ai) = regions{reg}; 
    ai = ai+1; 

    aovDat.encRet(ai) = 'bothReact';
    aovDat.count(ai)= sum(curReg==4);
    aovDat.reg(ai) = regions{reg}; 
    ai = ai+1; 

    aovDat.encRet(ai) = 'noReact'; 
    aovDat.count(ai)= sum(curReg==1);
    aovDat.reg(ai) = regions{reg}; 
    ai = ai+1; 
   

end

writetable(aovDat, ...
    ['G:\My Drive\GitHub\HpcAccConnectivityProject' ...
    '\reactivityPlotDat_encRetSplit.csv'])