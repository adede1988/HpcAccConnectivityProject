function [errorChans, metaDat] = getAllChanDat(chanFiles, sub)

metaDat = []; 

Ei = 1; 

errorChans = []; 

subIDs = cellfun(@(x) split(x, '_'), {chanFiles.name}, 'UniformOutput',false);
subIDs_all = cellfun(@(x) x{2}, subIDs, 'UniformOutput',false); 
subIDs_uni = unique(subIDs_all); 
clear subIDs

% for sub = 1:length(subIDs_uni)
%     sub
    Ei = 1; 
    EiSave = Ei; 
    try

    %% load data 
    curSub = chanFiles(cellfun(@(x) strcmp(x, subIDs_uni{sub}), subIDs_all)); 
    curSubReact = zeros(length(curSub), 1); %indicator variable for HFB reactive channels
    chanDat = load([curSub(1).folder '/' curSub(1).name]).chanDat; 
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

    if acc>0 && chanDat.age > 13 %memory and age filters

    if isfile(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' chanDat.site '_' chanDat.subID '_all' '.mat'])
        metaDat = load(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' chanDat.site '_' chanDat.subID '_all' '.mat']).metaDat;
%         metaDat.leadLag = rmfield(metaDat.leadLag, {'subMem', 'retMem'});
    else
%     curIdx = []; 
    flag = true;
    for chan = 1:length(curSub)
        chanDat = load([curSub(chan).folder '/' curSub(chan).name]).chanDat; 
        chanDat.reactive = 1; 
        chanDat.goodSub = 1; 
        chanDat.allBrod = 1; 
        if isfield(chanDat, 'HFB') && isfield(chanDat, 'leadLag4') && isfield(chanDat, 'ISPC')%check for complete processing
            %only reactive channels have leadlag calculations 
            %KNOWN ISSUE CORRECTED IN NEXT RUN: HFB for RT locked encoding
            %data had the wrong baseline. This caused spurious reactive
            %detections. inclChan in the leadLag3 struct was not
            %contaminated by this because HFB was recalculated on chan2
            %data (see singleChanPipeline in the lead Lag analysis block
            %Mostly channels that should not have been included had leadLag
            %calculations run. However, it did happen occasionally that a
            %channel which should have had calcualtions run was missed.
             if isstruct(chanDat.leadLag4)
                reactive = chanDat.leadLag4.inclChan; 
                if ismember(chan, reactive)
                    reactive = 1;
                else
                    reactive = 0; 
                end
             else
                reactive = 0; 
             end


             if sum(reactive>0)>0
                curSubReact(chan) = 1; 
%                 curIdx = [curIdx Ei]; 
                if Ei == 1
                    allDat = chanDat; 
                    Ei = Ei+1; 
                else
                    allDat(Ei) = chanDat; 
                    Ei = Ei + 1; 
                end

             end

        else
            %note down error channel! 
            flag = false; 
            errorChans = [errorChans sub*1000+chan];
            disp(['missing processing for: ' subIDs_uni{sub} ' channel: ' num2str(chan)])

        end

        

    end
    clear chanDat

    %% aggregating multichannel data
    if flag
        
        %key data for HFB, leadLag, TF, and ISPC
        %each one will be handled as a separate saved output
        %preparing to create new aggregate variables for all reactive
        %channels

        HFB = allDat(1).HFB;
        leadLag = allDat(1).leadLag4;  
        %hard code extra space for variable number of detected significant
        %clusters! 
%         leadLag.subMem(:,:,end+1:1000, :) = nan; 
%         leadLag.retMem(:,:,end+1:1000, :) = nan; 
        TF = allDat(1).TF; 
        ISPC = allDat(1).ISPC;
        HFB_names = fieldnames(HFB); 
        TF_names = fieldnames(TF); 
        ispc_names = fieldnames(ISPC); 
        leadLag_names = fieldnames(leadLag); 
        Nreact = sum(curSubReact); 
        %eliminate non-reactive channels from connectivity measures (ISPC and leadLag) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        for chan = 1:length(allDat)
            for fi = 1:length(leadLag_names)
                dimVals = size(allDat(chan).leadLag4.(leadLag_names{fi}));  
                if ismember(length(curSubReact), dimVals) %check that there's a channel dimension in the field
                    S.subs = repmat({':'}, 1, ndims(allDat(chan).leadLag4.(leadLag_names{fi})));
                    S.subs{1} = find(curSubReact==0);
                    S.type = '()';
                    allDat(chan).leadLag4.(leadLag_names{fi}) = subsasgn(allDat(chan).leadLag4.(leadLag_names{fi}), S, []); 
                end
            end
            for fi = 1:length(ispc_names)
                dimVals = size(allDat(chan).ISPC.(ispc_names{fi})); 
                if ismember(length(curSubReact), dimVals) %check that there's a channel dimension in the field
                    allDat(chan).ISPC.(ispc_names{fi})(curSubReact==0,:,:,:) = []; 
                end
            end
        end
    
        %expand the output variables to contain all channel data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %while also contracting the channel count to reactive channels
        for fi = 1:length(HFB_names)
            if ismember(fi, [1,2,5,6,9,10,11,12,15,16,17,18]) %hard code ugly for which HFB fields hold the output data
                HFB.(HFB_names{fi}) = zeros([sum(curSubReact), ...
                                        size(allDat(1).HFB.(HFB_names{fi}))]);
            end
        end

        for fi = 1:length(leadLag_names)
            dimVals = size(allDat(1).leadLag4.(leadLag_names{fi}));
            if ismember(sum(curSubReact), dimVals) %connectivity measures can avoid the hard code w/ channel count dimension
                leadLag.(leadLag_names{fi}) = zeros([sum(curSubReact), ...
                                                size(allDat(1).leadLag4.(leadLag_names{fi})) ]);
            end
        end

        for fi = 1:length(TF_names) %avoids field selection altogether because it only has data fields (best practice)
            TF.(TF_names{fi}) = zeros([sum(curSubReact), ...
                                        size(allDat(1).TF.(TF_names{fi}))]);     
        end

        for fi = 1:length(ispc_names)
            dimVals = size(allDat(1).ISPC.(ispc_names{fi}));
            if ismember(sum(curSubReact), dimVals)
                ISPC.(ispc_names{fi}) = zeros([sum(curSubReact), ...
                                            size(allDat(1).ISPC.(ispc_names{fi}) ) ] ); 
            end
        end

        
        %loop channels and store the values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              

        for chan = 1:length(allDat)
            %HFB
            for fi = 1:length(HFB_names)
                if ismember(fi, [1,2,5,6,9,10,11,12,15,16,17,18]) 
                    HFB.(HFB_names{fi})(chan,:,:) = allDat(chan).HFB.(HFB_names{fi});
                end
            end

            %leadLag
            for fi = 1:length(leadLag_names)
                dimVals = size(allDat(chan).leadLag4.(leadLag_names{fi}));
                if ismember(sum(curSubReact), dimVals) && fi>1 %skip the inclChan variable
                    tmp = allDat(chan).leadLag4.(leadLag_names{fi});
                    leadLag.(leadLag_names{fi})(chan,:,:,:,:) = tmp;
                end
            end

            %TF
            for fi = 1:length(TF_names) 
                TF.(TF_names{fi})(chan,:,:) = allDat(chan).TF.(TF_names{fi});     
            end






            %ISPC
            for fi = 1:length(ispc_names)
                dimVals = size(allDat(chan).ISPC.(ispc_names{fi}));
                if ismember(sum(curSubReact), dimVals)
                    ISPC.(ispc_names{fi})(chan,:,:,:,:) = allDat(chan).ISPC.(ispc_names{fi}); 
                end
            end
     
            %free up space
            allDat(chan).HFB = 1; 
            allDat(chan).leadLag3 = 1; 
            allDat(chan).leadLag2 = 1; 
            allDat(chan).leadLag4 = 1; 
            allDat(chan).ISPC = 1;
            allDat(chan).TF = 1; 
           
        end



  %connectivity was only calculated in one direction, so mirror the matrix
  %to copy the connectivity data (should be symmetric) 
for fi = 5:16
    ISPC.(ispc_names{fi}) = multiDmirror(ISPC.(ispc_names{fi}));
end


zachLabs = readtable('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject/brodmann_by_subj.csv');
zachLabs = zachLabs(cell2mat(cellfun(@(x) strcmp(x, allDat(1).subID), {zachLabs.subj}, 'uniformoutput', false )), :);
zachLabs(isnan(zachLabs.x), :) = []; 
zachLabs(curSubReact==0, :) = []; 
brod_new = table2array(zachLabs(:,5)); 
brod_old = {allDat.brodmann}; 

if ~strcmp(brod_old{1}, 'ERROR')
    disp(['match to old brod:' num2str(strcmp(brod_new{1}, brod_old{1})) ])
else
    disp(brod_new{1})
end
   
       
        metaDat = allDat(1);
        metanames = fieldnames(metaDat); 
        metaDat = rmfield(metaDat, metanames(19:end));
        metaDat.reactiveChans = curSubReact; 
        metaDat.elecpos = allDat(1).elecpos(curSubReact==1,:); 
        metaDat.meetLabs = allDat(1).labels(curSubReact==1,:);
        metaDat.brodmann = brod_new; 
        metaDat.chi = [allDat.chi];
        clear allDat

        metaDat.HFB = HFB; 
        metaDat.TF = TF; 
        metaDat.leadLag = leadLag; 
        metaDat.ISPC = ISPC; 

        %save the outputs one at a time: 
   
        save(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' metaDat.site '_' metaDat.subID '_all' '.mat'], 'metaDat', '-v7.3')
%         clear HFBout
% 
%         TFout = metaDat; 
%         TFout.TF = TF; 
%         clear TF
%         save(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' metaDat.site '_' metaDat.subID '_TF' '.mat'], 'TFout', '-v7.3')
%         clear TFout 
% 
%         ISPCout = metaDat; 
%         ISPCout.ISPC = ISPC; 
%         clear ISPC
%         save(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' metaDat.site '_' metaDat.subID '_ISPC' '.mat'], 'ISPCout', '-v7.3')
%         clear ISPCout 
% 
%         LLout = metaDat; 
%         LLout.LL = leadLag; 
%         clear leadLag
%         save(['R:\MSS\Johnson_Lab\dtf8829\SUMDAT\' metaDat.site '_' metaDat.subID '_LL' '.mat'], 'LLout', '-v7.3')
%         clear LLout 




    else
        subIDX = find(cellfun(@(x) strcmp(x, subIDs_uni{sub}), {allDat.subID})); 
        for chan = 1:length(subIDX)
            allDat(subIDX(chan)).goodSub = 0; 
        end
    end
    end
    end %memory and age filter end
    catch
        disp(['sub ' num2str(sub) ' ' curSub(chan).name 'fail'])
        metaDat = []; 
        Ei = EiSave; 

    end

   


end