function valMat = quickChanFillIn(valMat)

%HARD CODED function to add transpose to all connectivity matrices! 

allDims = size(valMat); 

for ii = 1:allDims(2)
    for jj = 1:allDims(3)
        for kk = 1:allDims(4)
            cur = squeeze(valMat(:,ii,jj,kk,:)); 
            cur = cur+ cur'; 
            valMat(:,ii,jj,kk,:) = cur; 
        end
    end
end





end

