function [padDat] = mirrorPad(dat)
        if size(dat,1) > size(dat,2)
        padDat = flip(dat,1);
        padDat = [padDat; dat; padDat];
        else
        padDat = flip(dat,2);
        padDat = [padDat, dat, padDat]; 
        end
end
