%audit all of the fig3Dat (connectivity) to find key information about all
%significant connections for further investigation

%figFile, encret, HFBimage, cc (cluster idx num), pp (channel pair idx num)



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

fig3Dat = dir([figDat 'Figure3']); 
fig3Dat(1:2) = []; 

phases = {'enc', 'ret'}; 
statName = {'HFB', 'image'};


parfor ii = 1:length(fig3Dat)
    
    curDat = load([fig3Dat(ii).folder '/' fig3Dat(ii).name]).outDat; 

    for er = 1:length(phases)
        for stat = 1:length(statName)
           
            fn = fig3Dat(ii).name; 
            connectionAudit(curDat, regions, phases{er}, ...
                statName{stat}, datPre, fn);

        
        end
    end
    


end