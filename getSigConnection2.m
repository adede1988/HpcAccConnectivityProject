function [] = getSigConnection(aggTargs, allDat)
% dimVals = size(allDat{3}.leadLag.subMem, [4,5]);
dimVals2 = size(allDat{3}.leadLag.retMem, [4,5]); 
%reg X reg X tVals / corr diff / p values X offset X time
% sigConSub = zeros([length(aggTargs), length(aggTargs), 3,  dimVals] );
% sigConRet = zeros([length(aggTargs), length(aggTargs), 3,  dimVals2] ); 
% path = 'G:\My Drive\Johnson\MTL_PFC_networkFigs\HFB_leadLag230910\';
for reg1 = 1:length(aggTargs)
    reg1
    for reg2 = 1:length(aggTargs)
        tic
        
        disp(['working on: ' num2str(reg1) ' X ' num2str(reg2)])
%         regRes = nan([2, dimVals, 100]);
        regRes2 = nan([2, dimVals2, 100]); %retrieval 
        regSubs = nan([100,1]);
        regSubIDs = cell(100,1); 
        chani = cell(100,1); 
        regd = regSubs; 
        regacc = regSubs; 
%         submissRT = regd; 
%         subhitRT = regd; 
        retmissRT = regd; 
        rethitRT = regd; 
        %at some point later it'd be better if these were done by trial
        %before averaging but that's not possible without redoing the
        %signal processing
%         regResLat = nan([2,100]); %latency of strongest connection (regardless of leadlag)
        regRes2Lat = nan([2,100]); 
%         regResOff = nan([2,dimVals(2), 100]); %calcualte weighted mean of the offset at each time point
        regRes2Off = nan([2,dimVals2(2), 100]); 
        ri = 1; 
        for sub = 1:length(allDat)
           

            
            if ~isempty(allDat{sub})


                T = sum(allDat{sub}.retInfo(:,1)==1 | allDat{sub}.retInfo(:,1)==2); 
                Hr = sum(allDat{sub}.retInfo(:,1)==1) / T; 
                T = sum(allDat{sub}.retInfo(:,1)==3 | allDat{sub}.retInfo(:,1)==4); 
                F = sum(allDat{sub}.retInfo(:,1)==4);
                if F == 0
                    acc = Hr - F/T; 
                    F = 1; 
                else
                    acc = Hr - F/T; 
                end
                Fr = F / T; 
               
                d = norminv(Hr) - norminv(Fr); 

%                 smRT = mean(allDat{sub}.encInfo(allDat{sub}.misses & allDat{sub}.use, 4));
%                 shRT = mean(allDat{sub}.encInfo(allDat{sub}.hits & allDat{sub}.use, 4));
                rmRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==2, 3));
                rhRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==1, 3));

                c = allDat{sub}.leadLag; 
                b = allDat{sub}.meetLabs(:,3); 
                ID = allDat{sub}.subID;
                reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
                reg2i = find(cellfun(@(x) sum(strcmp(aggTargs(reg2).lab, x)), b));

                if ~isempty(reg1i) && ~isempty(reg2i) 
                    for i1 = 1:length(reg1i)
                        for i2 = 1:length(reg2i)
                            if reg1i(i1) ~= reg2i(i2)
%                             regRes(:, :, :, ri) = c.subMem(reg1i(i1), reg2i(i2), :, :, :); 
                            regRes2(:,:, :, ri) = c.retMem(reg1i(i1), reg2i(i2), :, :, :); 
                            regSubs(ri) = sub; 
                            regSubIDs{ri} = ID; 
                            regd(ri) = d; 
                            regacc(ri) = acc; 
%                             submissRT(ri) = smRT; 
%                             subhitRT(ri) = shRT; 
                            retmissRT(ri) = rmRT; 
                            rethitRT(ri) = rhRT; 
                            chani{ri} = [num2str(reg1i(i1)) '_' num2str(reg2i(i2))];
                            
                            tim = allDat{sub}.leadLag.encTim; 
                            %sub Hit latency:                            
%                             regResLat(1,ri) = getLatency(squeeze(c.subMem(reg1i(i1), reg2i(i2), 1, :, :)), tim, shRT);  
                            %sub Miss latency: 
%                             regResLat(2,ri) = getLatency(squeeze(c.subMem(reg1i(i1), reg2i(i2), 2, :, :)), tim, smRT); 
                            tim = allDat{sub}.leadLag.retTim; 
                            %ret Hit latency: 
                            regRes2Lat(1,ri) = getLatency(squeeze(c.retMem(reg1i(i1), reg2i(i2), 1, :, :)), tim, rhRT);  
                            %sub Miss latency: 
                            regRes2Lat(2,ri) = getLatency(squeeze(c.subMem(reg1i(i1), reg2i(i2), 2, :, :)), tim, rmRT); 

                            
%                             regResOff(1,:,ri) = getOffset(squeeze(c.subMem(reg1i(i1), reg2i(i2), 1, :, :))); 
%                             regResOff(2,:,ri) = getOffset(squeeze(c.subMem(reg1i(i1), reg2i(i2), 2, :, :))); 
                            regRes2Off(1,:,ri) = getOffset(squeeze(c.retMem(reg1i(i1), reg2i(i2), 1, :, :))); 
                            regRes2Off(2,:,ri) = getOffset(squeeze(c.retMem(reg1i(i1), reg2i(i2), 2, :, :))); 
                            

          



                            ri = ri + 1; 
                            end
                        end
                    end


                end

               


            end
        end
        if ri < 100
%             regRes(:,:,:,ri:end) = []; 
            regRes2(:,:,:,ri:end) = [];
            regSubs(ri:end) = []; 
            regSubIDs(ri:end) = []; 
            regd(ri:end) = []; 
            regacc(ri:end) = []; 
            chani(ri:end) = [];
%             submissRT(ri:end) = []; 
%             subhitRT(ri:end) = []; 
            retmissRT(ri:end) = []; 
            rethitRT(ri:end) = [];
%             regResLat(:,ri:end) = []; %latency of strongest connection (regardless of leadlag)
            regRes2Lat(:,ri:end) = []; 
%             regResOff(:,:,ri:end) = []; %calcualte weighted mean of the offset at each time point
            regRes2Off(:,:,ri:end) = []; 
        end

        % package necessary data for cluster analysis
        
        LLdat = struct; 
%         LLdat.regRes = regRes; 
        LLdat.regRes2 = regRes2; 
        LLdat.regSubs = regSubs; 
        LLdat.regSubIDs = regSubIDs; 
        LLdat.aggTargs = aggTargs; 
        LLdat.reg1 = reg1; 
        LLdat.reg2 = reg2; 
        LLdat.d = regd; 
        LLdat.acc = regacc; 
        LLdat.n_sub = length(unique(regSubs));
        LLdat.n_pair = length(regSubs); 
        LLdat.encTim = allDat{3}.leadLag.encTim;  
        LLdat.retTim = allDat{3}.leadLag.retTim;  
%         LLdat.submissRT = submissRT; 
%         LLdat.subhitRT = subhitRT; 
        LLdat.retmissRT = retmissRT; 
        LLdat.rethitRT = rethitRT;
        LLdat.chani = chani; 
%         LLdat.regResLat = regResLat; %latency of strongest connection (regardless of leadlag)
        LLdat.regRes2Lat = regRes2Lat; 
%         LLdat.regResOff = regResOff; %calcualte weighted mean of the offset at each time point
        LLdat.regRes2Off = regRes2Off; 
            

        save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_LL_KEY_STATS\' aggTargs(reg1).lab '_' aggTargs(reg2).lab '_retrieve.mat'], 'LLdat', '-v7.3')


      



        disp(['........................' num2str(round(toc/60, 1))])
    end
end






end















