function [] = visualizeTFsingleTrialDat(allRes, reg, regions, phase, fi)


    test = cellfun(@(x) strcmp(x, regions{reg}), {allRes{:,5}});
    curReg = allRes(test,:);
    hits = [];
    misses = []; 
    hitRT = [];
    missRT = []; 
    hitLat = []; 
    missLat = []; 
    for ii = 1:size(curReg,1)
        hits = [hits, curReg{ii,9+(fi-1)*2}];
        hitRT = [hitRT; curReg{ii, 3}]; 
        misses = [misses, curReg{ii,10+(fi-1)*2}];
        missRT = [missRT; curReg{ii, 4}]; 
        hitLat = [hitLat; curReg{ii, 13}];
        missLat = [missLat; curReg{ii, 14}];

    end



    tim = curReg{1,6}; 

    hits= hits'; %reshape(hits, [flip(size(hits)), 1]);
    misses= misses'; %reshape(misses, [flip(size(misses)), 1]);

%eliminate bad trials
hits(hitLat==-1, :) = []; 
hitRT(hitLat==-1) = []; 
hitLat(hitLat==-1) = []; 

misses(missLat==-1, :) = []; 
missRT(missLat==-1) = []; 
missLat(missLat==-1) = [];

if fi==1
    freqLab = 'low'; 
else
    freqLab = 'high';
end


    quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
        phase, regions{reg}, 1, hitRT, missRT, freqLab)


end