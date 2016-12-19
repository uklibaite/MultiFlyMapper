function [imageOut,mask] = segmentImage_combo(image,alignInfo)

    
    s = size(image);
    midPoint = round(s./2);
    

    if ~isfield(alignInfo,'dilateSize') || isempty(alignInfo.dilateSize)
        dilateSize = 3;
    else
        dilateSize = alignInfo.dilateSize;
    end
    
    
    if ~isfield(alignInfo,'cannyParameter') || isempty(alignInfo.cannyParameter)
        cannyParameter = .1;
    else
        cannyParameter = alignInfo.cannyParameter;
    end
    
    
    if ~isfield(alignInfo,'imageThreshold') || isempty(alignInfo.imageThreshold)
        threshold = 3;
    else
        threshold = alignInfo.imageThreshold;
    end
    
    
    if ~isfield(alignInfo,'minimumArea') || isempty(alignInfo.minimumArea)
        minimumArea = 5000;
    else
        minimumArea = alignInfo.minimumArea;
    end
    
   
    if ~isfield(alignInfo,'maxDilateSize') || isempty(alignInfo.maxDilateSize)
        maxDilateSize = 6;
    else
        maxDilateSize = alignInfo.maxDilateSize;
    end


    minCannyParameter = 0;
    
    image2 = imcomplement(image);
    image2(image2 < threshold) = 0;
    
    E = edge(image2,'canny',cannyParameter,'nothinning');
    se = strel('square',dilateSize);
    E2 = imdilate(E,se);
    mask = imfill(E2,'holes');
    
    
    CC = bwconncomp(mask,8);
    if length(CC.PixelIdxList) > 1
        
        centroidVals = zeros(CC.NumObjects,1);
        for i=1:CC.NumObjects
            [ii,jj] = ind2sub(s,CC.PixelIdxList{i});
            ii = mean(ii);
            jj = mean(jj);
            centroidVals(i) = sum((midPoint - [ii jj]).^2);
        end
        
        idx = argmin(centroidVals);
        temp = false(size(mask));
        temp(CC.PixelIdxList{idx}) = true;
        mask = temp;
        
        %         idx = argmax(returnCellLengths(CC.PixelIdxList));
        %         temp = false(size(mask));
        %         temp(CC.PixelIdxList{idx}) = true;
        %         mask = temp;
    end
    
    
    while sum(mask(:)) < minimumArea && dilateSize <= maxDilateSize && cannyParameter > minCannyParameter
        
        dilateSize = dilateSize + 1;
        cannyParameter = .1;
        se = strel('square',dilateSize);
        E2 = imdilate(E,se);
        mask = imfill(E2,'holes');
        
        CC = bwconncomp(mask,8);
        if length(CC.PixelIdxList) > 1
            
            centroidVals = zeros(CC.NumObjects,1);
            for i=1:CC.NumObjects
                [ii,jj] = ind2sub(s,CC.PixelIdxList{i});
                ii = mean(ii);
                jj = mean(jj);
                centroidVals(i) = sum((midPoint - [ii jj]).^2);
            end
            
            idx = argmin(centroidVals);
            temp = false(size(mask));
            temp(CC.PixelIdxList{idx}) = true;
            mask = temp;
            
            %idx = argmax(returnCellLengths(CC.PixelIdxList));
            %temp = uint8(zeros(size(mask)));
            %temp(CC.PixelIdxList{idx}) = true;
            %mask = temp;
        end
        
        
    end
       
    
    if ~isinteger(image2)
        image2 = uint8(image2);
    end
    
    if ~islogical(mask)
        mask = mask > 0;
    end

    imageOut = immultiply(mask,image2);
    
    
    
    