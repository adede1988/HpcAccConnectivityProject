function [difVal] = myAngDif(ang1, ang2)

    difVal = ang1 - ang2;
    if abs(difVal) > pi
        difVal = (2*pi - abs(difVal)) * (-difVal/abs(difVal)); 
    end

   





end