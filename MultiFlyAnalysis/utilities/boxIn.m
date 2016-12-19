function [xs,ys] = boxIn(xs,ys,s)

    if xs(1) < 1
        xs = xs - xs(1) + 1;
    else
        if xs(end) > s(2)
            xs = xs + (s(2) - xs(end));
        end
    end
    
    
    if ys(1) < 1
        ys = ys - ys(1) + 1;
    else
        if ys(end) > s(1)
            ys = ys + (s(1) - ys(end));
        end
    end