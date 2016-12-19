function out = keepLargestConnectedComponent(image)

    CC = bwconncomp(image>0);
    
    if CC.NumObjects > 1
        
        idx = argmax(returnCellLengths(CC.PixelIdxList));
        mask = false(size(image));
        mask(CC.PixelIdxList{idx}) = true;
        
        out = image;
        out(~mask) = 0;
        
    else
        
        out = image;
        
    end