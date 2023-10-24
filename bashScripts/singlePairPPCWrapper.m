

%github access token: 
%  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
%  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
% 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
% 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


%simulated HPC looping input from shell script

%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';


%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])

savFolder = [datPre 'PPC/'];
%% initialize chanFiles, dif versions for local v. quest running


% chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;
chanFiles = load([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat']).chanFiles;

chanFiles(~[chanFiles.targPair]) = []; 


% % how the partnerChan files were made: 
% datFolder = [datPre 'CHANDAT']; 
% chanFiles = dir(datFolder);
% test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
% chanFiles = chanFiles(test); 


%I need a way of assigning the partner channel 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi <length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = [chanFiles(chan).folder '/CHANRAW']; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = subFiles(chan2).folder; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end
% 
% save([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat'], 'partnerChans')
% 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi < length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = '/projects/p31578/dtf8829/QuestConnect/CHANDAT/CHANRAW'; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = '/projects/p31578/dtf8829/QuestConnect/CHANDAT/CHANRAW'; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end
% 
% save([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat'], 'partnerChans')
% % HFBdat = load('R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\hip.mat').HFBdat; 
% % regions = {HFBdat.aggTargs.lab}; 
% % regions(2) = []; 
% % highfrex = linspace(70, 150, 81); 
% % tic
% % parfor ii = 1:length(chanFiles)
% %     ii
% % %     if isempty(chanFiles(ii).targPair)
% %     chanDat2 = load([chanFiles(ii).folder '/' chanFiles(ii).name2]).chanDat; 
% % 
% %     chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
% % 
% %     targ1 = sum(cellfun(@(x) strcmp(chanDat.labels(chanDat.chi, 3),x), regions)) > 0;
% %     targ2 = sum(cellfun(@(x) strcmp(chanDat2.labels(chanDat2.chi, 3),x), regions)) > 0;
% %     if targ1 && targ2
% %         chanFiles(ii).targPair = true; 
% %         chanFiles(ii).chan1Reg = chanDat.labels(chanDat.chi,3); 
% %         chanFiles(ii).chan2Reg = chanDat2.labels(chanDat2.chi,3); 
% % %         %both in target regions, so check HFB reactivity
% % %         HFB1 = getHFB(chanDat, highfrex); 
% % %         if sum(reactiveTest_100(HFB1)==1)>0
% % %             HFB2 = getHFB(chanDat2, highfrex);
% % %             if sum(reactiveTest_100(HFB2)==1)>0
% % %                 chanFiles(ii).bothReact = true; 
% % %             else
% % %                 chanFiles(ii).bothReact = false; 
% % %             end
% % %         else
% % %             chanFiles(ii).bothReact = false; 
% % %         end
% %     else
% %         chanFiles(ii).targPair = false; 
% %         chanFiles(ii).chan1Reg = chanDat.labels(chanDat.chi,3); 
% %         chanFiles(ii).chan2Reg = chanDat2.labels(chanDat2.chi,3); 
% %     end
% %         
% % %     end
% % end
% % toc
% % save([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat'], 'chanFiles')
% % 
% % tmp = chanFiles; 
% % chanFiles = load([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat']).partnerChans;
% % parfor ii = 1:length(chanFiles)
% %     
% %     chanFiles(ii).targPair = tmp(ii).targPair; 
% %     chanFiles(ii).chan1Reg = tmp(ii).chan1Reg; 
% %     chanFiles(ii).chan2Reg = tmp(ii).chan2Reg; 
% % 
% % 
% % end
% % save([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat'], 'chanFiles')

%% run the pipeline


disp(['going for ' chanFiles(start).name ' and ' chanFiles(start).name2] )
tic
singlePairPPCpipeline(chanFiles(start), savFolder); 
toc

