function [w] = basis(x, y, xc, yc, hc)

w = 0;
diffx = x - xc;
diffy = y - yc;

if abs(diffx) >= hc | abs(diffy) >= hc
    return
end

if diffx >= 0 & diffy >= 0
    if diffx + diffy > hc
        return
    else
        w = (-diffx - diffy) / hc + 1;
    end

elseif diffx < 0 & diffy >= 0
    if diffx + diffy > 0
        w = -diffy / hc + 1;
    else
        w = diffx / hc + 1;
    end

elseif diffx < 0 & diffy < 0
    if diffx + diffy < -hc
        return
    else
        w = (diffx + diffy) / hc + 1;
    end

else
    if diffx + diffy > 0
        w = -diffx / hc + 1;
    else
        w = diffy / hc + 1;
    end
end