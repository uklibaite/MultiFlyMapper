function out = calculateKLCost(x,ydata,ps)


    d = findListDistances(x,ydata);
    d = d.^2;
    out = log(sum((1+d).^-1)) + sum(ps.*log(1+d));