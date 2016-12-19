function x = boundbox(x,low,high)

    if length(x) == 1

    if x < low
        x = low;
    else
        if x > high
            x = high;
        end
    end
    
    else
       
        x(x < low) = low;
        x(x > high) = high;
                
    end