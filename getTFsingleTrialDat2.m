function latency = getTFsingleTrialDat2(allRes, reg, regions, phase)
    statInfo = struct; %output for stats
    statInfo.reg = regions{reg}; 
    statInfo.regi = reg; 
    statInfo.phase = phase; 

    %create vector to identify the data for the target region 
    test = cellfun(@(x) strcmp(x, regions{reg}), {allRes{:,5}});
    curReg = allRes(test,:);
    hits = [];
    misses = []; 
    hits_p = []; 
    misses_p = []; 
    hitRT = [];
    missRT = []; 
    hitChi = []; 
    missChi = []; 
    hitTriali = []; 
    missTriali = []; 
    hitSub = cell(5000, 1); 
    missSub = cell(5000,1); 
    hitAge = []; 
    hitAcc = []; 
    hitd = []; 
    missAge = [];
    missAcc = [];
    missd = [];
    latency = []; 
    latency2 = []; 
    hi = 1; 
    mi = 1; 


    for ii = 1:size(curReg,1)
        hits = [hits, curReg{ii,9}];
        hits_p = [hits, curReg{ii,11}];
        hitRT = [hitRT; curReg{ii, 3}]; 
        misses = [misses, curReg{ii,10}];
        misses_p = [misses, curReg{ii,12}];
        missRT = [missRT; curReg{ii, 4}]; 
        latency = [latency; curReg{ii, 13}]; 
        latency2 = [latency2; curReg{ii, 14}]; 

        L = size(curReg{ii,9},2);
        hitChi = [hitChi; repmat(curReg{ii, 8}, L, 1) ]; 
        hitAge = [hitAge; repmat(curReg{ii, 21}, L, 1) ]; 
        hitAcc = [hitAcc; repmat(curReg{ii, 20}, L, 1) ]; 
        hitd = [hitd; repmat(curReg{ii, 19}, L, 1) ]; 
        hitTriali = [hitTriali; [1:L]' ]; 
        for jj = 0:L-1
            hitSub{hi+jj} = curReg{ii, 7}; 
        end
        hi = hi + L;
        L1 = L; 

        L = size(curReg{ii,10},2);
        missChi = [missChi; repmat(curReg{ii, 8}, L, 1) ];
        missAge = [missAge; repmat(curReg{ii, 21}, L, 1) ];
        missAcc = [missAcc; repmat(curReg{ii, 20}, L, 1) ];
        missd = [missd; repmat(curReg{ii, 19}, L, 1) ];
        missTriali = [missTriali; [L1+1:L1+L]' ];
        for jj = 0:L-1
            missSub{mi+jj} = curReg{ii, 7}; 
        end
        mi = mi + L; 
%         missSub = [missSub; repmat(curReg{ii, 7}, size(curReg{ii,2},2), 1) ];
    end
    hitSub(hi:end) = []; 
    missSub(mi:end) = []; 

    tim = curReg{1,6}; 
   

 
    eliminate = latency == -1;
    eliminate2 = latency2 == -1; 

    latency(eliminate) = []; 
    latency2(eliminate2) = []; 
    hits(:,eliminate, :) = [];
    hits_p(:,eliminate,:) = [];
    hitRT(eliminate) = []; 
    misses(:,eliminate2, :) = [];
    misses_p(:,eliminate2,:) = [];
    missRT(eliminate2) = []; 
    hitChi(eliminate) = []; 
    missChi(eliminate2) = []; 
    hitSub(eliminate) = []; 
    missSub(eliminate2) = []; 
    hitTriali(eliminate) = []; 
    missTriali(eliminate2) = []; 
    hitAge(eliminate) = []; 
    hitAcc(eliminate) = []; 
    hitd(eliminate) = []; 
    missAge(eliminate2) = [];
    missAcc(eliminate2) = [];
    missd(eliminate2) = [];



    statInfo.hits = hits; 
    statInfo.misses = misses; 
    statInfo.hits_p = hits_p;
    statInfo.misses_p = misses_p;
    statInfo.tim = tim; 
    statInfo.hitRT = hitRT; 
    statInfo.missRT = missRT;
    statInfo.hitLat = latency; 
    statInfo.missLat = latency2; 
    statInfo.hitChi = hitChi; 
    statInfo.missChi = missChi; 
    statInfo.hitSub = hitSub; 
    statInfo.missSub = missSub; 
    statInfo.hitTriali = hitTriali; 
    statInfo.missTriali = missTriali; 
    statInfo.hitAge = hitAge; 
    statInfo.hitAcc = hitAcc ;
    statInfo.hitd = hitd ;
    statInfo.missAge = missAge;
    statInfo.missAcc = missAcc;
    statInfo.missd = missd;
    statInfo.frex = logspace(log10(2),log10(80),100);

    save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\TF_singleTrial\' ...
        statInfo.reg '_' statInfo.phase '_all' '.mat'], "statInfo", '-v7.3')
    

    %% scratch below here 
    % originally used when there was a split happening by frequency

%     statInfo = struct; %output for stats
%     statInfo.reg = regions{reg}; 
%     statInfo.regi = reg; 
%     statInfo.phase = phase; 
% 
%     %create vector to identify the data for the target region 
%     test = cellfun(@(x) strcmp(x, regions{reg}), {allRes{:,5}});
%     curReg = allRes(test,:);
%     hits = [];
%     misses = []; 
%     hitRT = [];
%     missRT = []; 
%     hitChi = []; 
%     missChi = []; 
%     hitTriali = []; 
%     missTriali = []; 
%     hitSub = cell(5000, 1); 
%     missSub = cell(5000,1); 
%     hitAge = []; 
%     hitAcc = []; 
%     hitd = []; 
%     missAge = [];
%     missAcc = [];
%     missd = [];
%     latency = []; 
%     latency2 = []; 
%     hi = 1; 
%     mi = 1; 
% 
% 
%     for ii = 1:size(curReg,1)
%         hits = [hits, curReg{ii,11}];
%         hitRT = [hitRT; curReg{ii, 3}]; 
%         misses = [misses, curReg{ii,12}];
%         missRT = [missRT; curReg{ii, 4}]; 
%         latency = [latency; curReg{ii, 13}]; 
%         latency2 = [latency2; curReg{ii, 14}]; 
% 
%         L = size(curReg{ii,9},2);
%         hitChi = [hitChi; repmat(curReg{ii, 8}, L, 1) ]; 
%         hitAge = [hitAge; repmat(curReg{ii, 21}, L, 1) ]; 
%         hitAcc = [hitAcc; repmat(curReg{ii, 20}, L, 1) ]; 
%         hitd = [hitd; repmat(curReg{ii, 19}, L, 1) ]; 
%         hitTriali = [hitTriali; [1:L]' ]; 
%         for jj = 0:L-1
%             hitSub{hi+jj} = curReg{ii, 7}; 
%         end
%         hi = hi + L;
%         L1 = L; 
% 
%         L = size(curReg{ii,10},2);
%         missChi = [missChi; repmat(curReg{ii, 8}, L, 1) ];
%         missAge = [missAge; repmat(curReg{ii, 21}, L, 1) ];
%         missAcc = [missAcc; repmat(curReg{ii, 20}, L, 1) ];
%         missd = [missd; repmat(curReg{ii, 19}, L, 1) ];
%         missTriali = [missTriali; [L1+1:L1+L]' ];
%         for jj = 0:L-1
%             missSub{mi+jj} = curReg{ii, 7}; 
%         end
%         mi = mi + L; 
% %         missSub = [missSub; repmat(curReg{ii, 7}, size(curReg{ii,2},2), 1) ];
%     end
%     hitSub(hi:end) = []; 
%     missSub(mi:end) = []; 
% 
%     tim = curReg{1,6}; 
%    
% 
%  
%     eliminate = latency == -1;
%     eliminate2 = latency2 == -1; 
% 
%     latency(eliminate) = []; 
%     latency2(eliminate2) = []; 
%     hits(:,eliminate) = [];
%     hitRT(eliminate) = []; 
%     misses(:,eliminate2) = []; 
%     missRT(eliminate2) = []; 
%     hitChi(eliminate) = []; 
%     missChi(eliminate2) = []; 
%     hitSub(eliminate) = []; 
%     missSub(eliminate2) = []; 
%     hitTriali(eliminate) = []; 
%     missTriali(eliminate2) = []; 
%     hitAge(eliminate) = []; 
%     hitAcc(eliminate) = []; 
%     hitd(eliminate) = []; 
%     missAge(eliminate2) = [];
%     missAcc(eliminate2) = [];
%     missd(eliminate2) = [];
% 
% 
% 
%     statInfo.hits = hits; 
%     statInfo.misses = misses; 
%     statInfo.tim = tim; 
%     statInfo.hitRT = hitRT; 
%     statInfo.missRT = missRT;
%     statInfo.hitLat = latency; 
%     statInfo.missLat = latency2; 
%     statInfo.hitChi = hitChi; 
%     statInfo.missChi = missChi; 
%     statInfo.hitSub = hitSub; 
%     statInfo.missSub = missSub; 
%     statInfo.hitTriali = hitTriali; 
%     statInfo.missTriali = missTriali; 
%     statInfo.hitAge = hitAge; 
%     statInfo.hitAcc = hitAcc ;
%     statInfo.hitd = hitd ;
%     statInfo.missAge = missAge;
%     statInfo.missAcc = missAcc;
%     statInfo.missd = missd;
% 
%     save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\TF_singleTrial\' ...
%         statInfo.reg '_' statInfo.phase '_high' '.mat'], "statInfo")
% 











end