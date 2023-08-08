function out = getExampleConnection(targBrod, allDat, fromReg, targTim, targOff, targTim2, targOff2)

targBrod(cellfun(@(x) strcmp('ERROR', x), {targBrod.lab})) = []; 
%area X area X time X offset X up/down
out = 5; 
tim = allDat{1}.leadLag.encTim; 
LLtim = -150:150; 
for sub = 1:length(allDat)
    if ~isempty(allDat{sub})
        c = allDat{sub}.leadLag; 
        b = allDat{sub}.brodmann; 
        for chan1 = 1:size(c.subMem,1)
            chan1i = cellfun(@(x) strcmp(b{chan1}, x), {targBrod.lab});
            if ismember(find(chan1i), fromReg) %this is one of the target channels! 
                
            for chan2 = 1:size(c.subMem,1)
                if chan2 ~= chan1
                chan2i = cellfun(@(x) strcmp(b{chan2}, x), {targBrod.lab});
                if sum(chan2i)>0 %this is one of the target channels! 
                    
 
                    for cc = 1:size(c.subMem, 4)
                        if ~isnan(c.subMem(chan1, chan2, 1, cc, 1)) %&& length(find(tim>=c.subMem(chan1, chan2,1,cc,3) & tim<= c.subMem(chan1,chan2,1,cc,4)))>1
                            if c.subMem(chan1, chan2, 1, cc, 1) > 0
                                
                               
                                if targTim>=c.subMem(chan1, chan2,1,cc,3) && targTim<=c.subMem(chan1,chan2,1,cc,4)
                                    if targOff>=c.subMem(chan1, chan2,1,cc,6) && targOff<=c.subMem(chan1, chan2,1,cc,7)

                                chi_1 = getChiString(allDat{sub}.chi(chan1));
                                chi_2 = getChiString(allDat{sub}.chi(chan2));


                                chan1Dat = load(['R:\MSS\Johnson_Lab\dtf8829\CHANDAT' '\chanDat_' allDat{sub}.subID '_' chi_1 '.mat']).chanDat;
                                chan2Dat = load(['R:\MSS\Johnson_Lab\dtf8829\CHANDAT' '\chanDat_' allDat{sub}.subID '_' chi_2 '.mat']).chanDat;

                                highfrex = linspace(70, 150, 81); 
                                    %CHAN 1
                                [pow, mulTim, mulFrex] = getChanMultiTF(chan1Dat.enc, highfrex, chan1Dat.fsample, chan1Dat.enctim, 1);  
                                pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
                                pow = cell2mat(pow); %organize
                                highnumfrex = length(mulFrex); 
                                pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
                                pow = squeeze(mean(pow, 3)); %take the mean over frequencies

                                HFB1 = pow; 
                                clear pow

                                    %CHAN 2
                                [pow, mulTim, mulFrex] = getChanMultiTF(chan2Dat.enc, highfrex, chan2Dat.fsample, chan2Dat.enctim, 1);  
                                pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
                                pow = cell2mat(pow); %organize
                                highnumfrex = length(mulFrex); 
                                pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
                                pow = squeeze(mean(pow, 3)); %take the mean over frequencies
                                pow2 = pow; 
                                clear pow

                                missidx = find(chan1Dat.use & chan1Dat.misses); 
                                hitidx = find(chan1Dat.use & chan1Dat.hits); 

                                dstim = chan1Dat.retOtim(501:25:end-500); 
                                alltim = chan1Dat.enctim; 

                                missTemp = zeros(length(missidx), 301, length(dstim)); 
                                hitTemp = zeros(length(hitidx), 301, length(dstim)); 
                                
                             
                                for offSet = -150:150 %negative means current Channel leads, positive means other channel leads
                                    
                                    if offSet<0
                                        HFB2 = [pow2(abs(offSet)+1:end,:); zeros([abs(offSet), size(pow2,2)] ) ] ;
                                    elseif offSet>0
                                        HFB2 = [zeros([abs(offSet), size(pow2,2)] );  pow2(1:end-abs(offSet),:)];
                                    end
                                    %sub miss 
                                 
                        
                                    missTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, missidx), HFB2(x-500:x+500, missidx)), ...
                                                                       dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(missidx), length(dstim)] );
                        
                        
                                    hitTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, hitidx), HFB2(x-500:x+500, hitidx)), ...
                                                                       dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(hitidx), length(dstim)] );
                        
                                 
                            
                                end



                                figure('position', [0, 0, 600, 1000], 'visible', false)
                                subplot 311
                                imagesc(dstim, LLtim, squeeze(mean(hitTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim, targOff, 50, 'red', 'filled')
                                title('subHit')
                                colorbar
                                subplot 312
                                imagesc(dstim, LLtim, squeeze(mean(missTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim, targOff, 50, 'red', 'filled')
                                title('subMiss')
                                colorbar
                                subplot 313
                                imagesc(dstim, LLtim, squeeze(mean(hitTemp,1))-squeeze(mean(missTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim, targOff, 50, 'red', 'filled')
                                title([chan1Dat.subID ' ' chan1Dat.brodmann ' to ' chan2Dat.brodmann])
                                colorbar
                                set(gcf, 'color', 'w')
                                export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\exampleConnections\' ...
                                    chan1Dat.subID '_' num2str(chan1Dat.chi) '_' num2str(chan2Dat.chi) '.jpg'])



                                    end
                                end                         
                            end
    


                        end

                        if ~isnan(c.subMem(chan1, chan2, 2, cc, 2)) %&& length(find(tim>=c.subMem(chan1, chan2,2,cc,3) & tim<= c.subMem(chan1,chan2,2,cc,4)))>1
                           if c.subMem(chan1, chan2, 2, cc, 1) > 0
                           

                               if targTim2>=c.subMem(chan1, chan2,2,cc,3) && targTim2<=c.subMem(chan1,chan2,2,cc,4)
                                    if targOff2>=c.subMem(chan1, chan2,2,cc,6) && targOff2<=c.subMem(chan1, chan2,2,cc,7)

                                chi_1 = getChiString(allDat{sub}.chi(chan1));
                                chi_2 = getChiString(allDat{sub}.chi(chan2));


                                chan1Dat = load(['R:\MSS\Johnson_Lab\dtf8829\CHANDAT' '\chanDat_' allDat{sub}.subID '_' chi_1 '.mat']).chanDat;
                                chan2Dat = load(['R:\MSS\Johnson_Lab\dtf8829\CHANDAT' '\chanDat_' allDat{sub}.subID '_' chi_2 '.mat']).chanDat;

                                highfrex = linspace(70, 150, 81); 
                                    %CHAN 1
                                [pow, mulTim, mulFrex] = getChanMultiTF(chan1Dat.enc, highfrex, chan1Dat.fsample, chan1Dat.enctim, 1);  
                                pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
                                pow = cell2mat(pow); %organize
                                highnumfrex = length(mulFrex); 
                                pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
                                pow = squeeze(mean(pow, 3)); %take the mean over frequencies

                                HFB1 = pow; 
                                clear pow

                                    %CHAN 2
                                [pow, mulTim, mulFrex] = getChanMultiTF(chan2Dat.enc, highfrex, chan2Dat.fsample, chan2Dat.enctim, 1);  
                                pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
                                pow = cell2mat(pow); %organize
                                highnumfrex = length(mulFrex); 
                                pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
                                pow = squeeze(mean(pow, 3)); %take the mean over frequencies
                                pow2 = pow; 
                                clear pow

                                missidx = find(chan1Dat.use & chan1Dat.misses); 
                                hitidx = find(chan1Dat.use & chan1Dat.hits); 

                                dstim = chan1Dat.retOtim(501:25:end-500); 
                                alltim = chan1Dat.enctim; 

                                missTemp = zeros(length(missidx), 301, length(dstim)); 
                                hitTemp = zeros(length(hitidx), 301, length(dstim)); 
                                
                             
                                for offSet = -150:150 %negative means current Channel leads, positive means other channel leads
                                    
                                    if offSet<0
                                        HFB2 = [pow2(abs(offSet)+1:end,:); zeros([abs(offSet), size(pow2,2)] ) ] ;
                                    elseif offSet>0
                                        HFB2 = [zeros([abs(offSet), size(pow2,2)] );  pow2(1:end-abs(offSet),:)];
                                    end
                                    %sub miss 
                                 
                        
                                    missTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, missidx), HFB2(x-500:x+500, missidx)), ...
                                                                       dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(missidx), length(dstim)] );
                        
                        
                                    hitTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, hitidx), HFB2(x-500:x+500, hitidx)), ...
                                                                       dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(hitidx), length(dstim)] );
                        
                                 
                            
                                end



                                figure('position', [0, 0, 600, 1000], 'visible', false)
                                subplot 311
                                imagesc(dstim, LLtim, squeeze(mean(hitTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim2, targOff2, 50, 'red', 'filled')
                                title('subHit')
                                colorbar
                                subplot 312
                                imagesc(dstim, LLtim, squeeze(mean(missTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim2, targOff2, 50, 'red', 'filled')
                                title('subMiss')
                                colorbar
                                subplot 313
                                imagesc(dstim, LLtim, squeeze(mean(hitTemp,1))-squeeze(mean(missTemp,1)))
                                yticks([-150:50:150])
                                yticklabels([150,100,50,0,50,100,150])
                                ylabel([chan2Dat.brodmann ' leads           ' chan1Dat.brodmann ' leads'])
                                hold on 
                                scatter(targTim2, targOff2, 50, 'red', 'filled')
                                title([chan1Dat.subID ' ' chan1Dat.brodmann ' to ' chan2Dat.brodmann])
                                colorbar
                                set(gcf, 'color', 'w')
                                export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\exampleConnections\zzMissGreater_' ...
                                    chan1Dat.subID '_' num2str(chan1Dat.chi) '_' num2str(chan2Dat.chi) '.jpg'])



                                    end
                                end    

                           end
    
                        end
                    end
                end
                end
            end
            end
        end
    




    end


end












end