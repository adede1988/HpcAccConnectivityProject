


%looking at lead Lag results
% load('R:\MSS\Johnson_Lab\dtf8829\HFBCONDAT\HFBCONDAT_IR84.mat')

addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\myFrequentUse\export_fig_repo')

LLfiles = dir('R:\MSS\Johnson_Lab\dtf8829\HFBCONDAT');
LLfiles(1:2) = []; 

%chan.chan  X time X offset (averaged across trials for each connection)
allConnections = struct; 
allConnections.dlpfc_hip = zeros(2000, 901, 41); 
allConnections.dlpfc_phg = zeros(2000, 901, 41); 
allConnections.dlpfc_acc = zeros(2000, 901, 41); 
allConnections.hip_phg = zeros(2000, 901, 41); 
allConnections.hip_acc = zeros(2000, 901, 41); 
allConnections.phg_acc = zeros(2000, 901, 41); 

conNames = fieldnames(allConnections); 

conCodes = [[1,2]; [1,3]; [1,4]; [2,3]; [2,4]; [3,4]];


for sub = 1:length(LLfiles)
    sub
     %subject failed in HPC because of memory
    curSub = load([LLfiles(sub).folder '/' LLfiles(sub).name]).curSub; 
    if isfield(curSub, "leadLagHit")
    for con = 1:length(conNames)
        reg1 = find(curSub(1).activeRoiMat(:,conCodes(con,1))==1); 
        reg2 = find(curSub(1).activeRoiMat(:,conCodes(con,2))==1); 
        if ~isempty(reg1) && ~isempty(reg2)
            
            for ii = 1:length(reg1)
                for jj = 1:length(reg2)
                    chani = find(allConnections.(conNames{con})(:,500,20) == 0 ,1 );
                    allConnections.(conNames{con})(chani, :, :) = squeeze(mean(curSub(1).leadLagHit(reg1(ii), reg2(jj), :, :, :), 3)); 
                    chani
                end
            end






        end





    end

    end
end

%% plotting 

for con = 1:length(conNames)
    chani = find(allConnections.(conNames{con})(:,500,20) == 0 ,1 );
    for ii = 1:chani-1
    figure('visible', 'off', 'position', [1,1,800,600])
    test = squeeze(allConnections.(conNames{con})(ii, :, :)); 
    reg1 = split(conNames{con}, '_');
    reg2 = reg1{2}; 
    reg1 = reg1{1}; 
    imagesc(test')
    yticks([1:5:41])
    yticklabels([-100:25:100])
    ylabel(['<-' reg1 ' leads                     ' reg2 ' leads->'])
    xticks([100:100:900])
    xticklabels(curSub(1).HFB.encMulTim([101:100:901]))
    yline(21, '--', 'color', 'red', 'linewidth', 3)
    xline(201, '--', 'color', 'green', 'linewidth', 3)
    colorbar
    caxis([-.3, .3])
    title([conNames{con} ' connection:' num2str(ii)], 'interpreter', 'none')
    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\HFBleadLag\' conNames{con} num2str(ii) '.jpg'])
    end

end


for con = 1:length(conNames)
    figure('position', [1,1,800,600])
    test = squeeze(mean(allConnections.(conNames{con}), 1)); 
    imagesc(test')
    reg1 = split(conNames{con}, '_');
    reg2 = reg1{2}; 
    reg1 = reg1{1}; 
    imagesc(test')
    yticks([1:5:41])
    yticklabels([-100:25:100])
    ylabel(['<-' reg1 ' leads                     ' reg2 ' leads->'])
    xticks([100:100:900])
    xticklabels(curSub(1).HFB.encMulTim([101:100:901]))
    yline(21, '--', 'color', 'red', 'linewidth', 3)
    xline(201, '--', 'color', 'green', 'linewidth', 3)
    colorbar
    caxis([-.01, .01])
    title([conNames{con} ' mean lead lag'], 'interpreter', 'none')
    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\HFBleadLag\' conNames{con} 'MEAN' '.jpg'])

end



%% scratch below here

    %dlpfc_hip
    reg1 = find(curSub(1).activeROIMat(:,))



    hipChans = find(curSub(1).activeRoiMat(:,2)==1);
    phgChans = find(curSub(1).activeRoiMat(:,3)==1);

    %check that both chans are represented
    if ~isempty(hipChans) && ~isempty(phgChans)

        for ii = 1:length(hipChans)
            for jj = 1:length(phgChans)

                hip_phg(chani, :, :) = squeeze(mean(curSub(1).leadLagHit(hipChans(ii), phgChans(jj), :, :, :), 3)); 
                chani = chani + 1; 
            end
        end



    end
    


    end
end

chani
for ii = 1:chani-1
    figure
    test = squeeze(hip_phg(ii,ii, :, :)); 
    imagesc(test')
    yticks([1:5:41])
    yticklabels([-100:25:100])
    ylabel('<-phg leads                     hip leads->')
    xticks([100:100:900])
    xticklabels(curSub(1).HFB.encMulTim([101:100:901]))
    yline(21, '--', 'color', 'red', 'linewidth', 3)
    xline(201, '--', 'color', 'green', 'linewidth', 3)
    colorbar

end

figure
test = squeeze(mean(hip_phg(1:chani-1, 1:chani-1, :, :), [1,2]));



