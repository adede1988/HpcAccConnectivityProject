function [rLabs] = roiLabel(rawLabs)

%goal is to select a set of ROI labels that encompass potentially multiple
%different input labels

dlPFC = {'sfs', 'sfg', 'mfg'};
acc = {'acc'}; 
hip = {'hip'}; 
phg = {'phg'}; 


rLabs = cell(size(rawLabs)); 
for li = 1:length(rawLabs)
    ch = rawLabs{li}; 

    if sum(cellfun(@(x) strcmp(ch, x), dlPFC)) > 0
                rLabs{li} = 'dlPFC'; 
    elseif sum(cellfun(@(x) strcmp(ch, x), acc)) > 0
                rLabs{li} = 'acc'; 
    elseif sum(cellfun(@(x) strcmp(ch, x), hip)) > 0
                rLabs{li} = 'hip'; 
    elseif sum(cellfun(@(x) strcmp(ch, x), phg)) > 0
                rLabs{li} = 'phg'; 
    else
        rLabs{li} = 'ZZZ'; 
    end


end














end