% connectivity integrate script

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])





%% TF final evaluation 
%it's assumed here that the TF quest pipeline has already been run and that
%the outputs are available 
%stat: 0 = HFB 1 = image
%TF power files
outStatFiles = dir([datPre 'pairFiles\out']);
outStatFiles(1:2, :) = []; 
for ii = 1:length(outStatFiles)
    splitName = split(outStatFiles(ii).name, '_'); 
    outStatFiles(ii).reg1 = splitName{1}; 
    outStatFiles(ii).reg2 = splitName{2}; 
    outStatFiles(ii).phase = splitName{3}; 
    splitName = split(splitName{4}, 'stat'); 
    outStatFiles(ii).stat = splitName{2}; 


end


%loop reg X reg doing aggregated stats
parfor regii = 1:9
    for regjj = 1:9
        disp(['regii:' num2str(regii) ' regjj:' num2str(regjj)])
        try
%         connectionIntegrate(outStatFiles, regii, regjj, 'ret', 0)
        connectionIntegrate(outStatFiles, regii, regjj)
        catch
        end
    end
end