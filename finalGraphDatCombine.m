
%figures should not be docked: 
set(0, 'defaultfigurewindowstyle', 'normal')
%local paths: 


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';
figDat = 'R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\';

% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

regions = {'acc', 'dlPFC', 'hip', ...
    'lTemp', 'iTemp', 'mtl', 'pcc', 'pPFC', 'vis'}; 


%graph based connectivity metrics getting a .csv for R analysis

graphDat = dir([datPre 'graphAnalysis/out']);
L = length(graphDat); 
graphDat(1:2) = [];
ii = 1;
cur = load([graphDat(ii).folder '/' graphDat(ii).name]).outDat;
graphDat(ii).reg1 = cur.reg1; 
graphDat(ii).reg2 = cur.reg2; 
graphDat(ii).timeSet = cur.HFB_Image;
graphDat(ii).time = cur.time; 
graphDat(ii).freq = cur.freq; 
graphDat(ii).encRet = cur.encRet; 
graphDat(ii).hitBC = cur.hitBC(cur.chi); 
graphDat(ii).missBC = cur.missBC(cur.chi); 
graphDat(ii).hitST = cur.hitST(cur.chi); 
graphDat(ii).missST = cur.missST(cur.chi); 
for ii = 1:L
    tic
    cur = load([graphDat(ii).folder '/' graphDat(ii).name]).outDat;

    cur.hitSparce = cur.hitMat; 
    cur.missSparce = cur.missMat; 
    cur.hitSparce(cur.hitMat<.1) = 0; 
    cur.missSparce(cur.missMat<.1) = 0; 
    cur.hitSparce(cur.hitSparce>0) = 1; 
    cur.missSparce(cur.missSparce>0) = 1; 

    hitPos = cur.hitMat; 
    hitPos(hitPos<0) = 0; 
    missPos = cur.missMat; 
    missPos(missPos<0) = 0; 

    
    graphDat(ii).reg1 = cur.reg1; 
    graphDat(ii).reg2 = cur.reg2; 
    graphDat(ii).subID = cur.subID; 
    graphDat(ii).chi = cur.chi; 
    graphDat(ii).chi2 = cur.chi2; 
    graphDat(ii).timeSet = cur.HFB_Image;
    graphDat(ii).time = cur.time; 
    graphDat(ii).freq = cur.freq; 
    graphDat(ii).encRet = cur.encRet; 
    graphDat(ii).hitBC = cur.hitBC(cur.chi); 
    graphDat(ii).missBC = cur.missBC(cur.chi); 
    graphDat(ii).hitST = cur.hitST(cur.chi); 
    graphDat(ii).missST = cur.missST(cur.chi); 
    graphDat(ii).hitSaturation = sum(cur.hitMat(cur.chi,:)>.1) / ...
                                    size(cur.hitMat,1); 
    graphDat(ii).missSaturation = sum(cur.missMat(cur.chi,:)>.1) / ...
                                    size(cur.hitMat,1);
    graphDat(ii).hitEff = efficiency_bin(cur.hitSparce); 
    graphDat(ii).missEff = efficiency_bin(cur.missSparce); 
    graphDat(ii).hitChar = charpath(distance_wei(1./hitPos)); 
    graphDat(ii).missChar = charpath(distance_wei(1./missPos)); 

    if strcmp(cur.HFB_Image, 'image') %also put the flip in there for image connections
        graphDat(ii+L).reg1 = cur.reg2; 
        graphDat(ii+L).reg2 = cur.reg1; 
        graphDat(ii+L).subID = cur.subID; 
        graphDat(ii+L).chi = cur.chi; 
        graphDat(ii+L).chi2 = cur.chi2; 
        graphDat(ii+L).timeSet = cur.HFB_Image;
        graphDat(ii+L).time = cur.time; 
        graphDat(ii+L).freq = cur.freq; 
        graphDat(ii+L).encRet = cur.encRet; 
        graphDat(ii+L).hitBC = cur.hitBC(cur.chi); 
        graphDat(ii+L).missBC = cur.missBC(cur.chi); 
        graphDat(ii+L).hitST = cur.hitST(cur.chi); 
        graphDat(ii+L).missST = cur.missST(cur.chi); 
        graphDat(ii+L).hitSaturation = sum(cur.hitMat(cur.chi,:)>.1) / ...
                                        size(cur.hitMat,1); 
        graphDat(ii+L).missSaturation = sum(cur.missMat(cur.chi,:)>.1) / ...
                                        size(cur.hitMat,1);
        graphDat(ii+L).hitEff = efficiency_bin(cur.hitSparce); 
        graphDat(ii+L).missEff = efficiency_bin(cur.missSparce); 
        graphDat(ii+L).hitChar = charpath(distance_wei(1./hitPos)); 
        graphDat(ii+L).missChar = charpath(distance_wei(1./missPos)); 


    end

toc
end

test = cellfun(@(x) isempty(x), {graphDat.reg1});
graphDat(test) = []; 

t = struct2table(graphDat);

% Save the table as a CSV file
writetable(t, [figDat 'graphDat.csv']);