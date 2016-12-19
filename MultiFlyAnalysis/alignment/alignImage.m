function [X,Y,rotationAngle,loopImage,asymValue] = alignImage(image,backgroundImage,alignInfo,angleGuess,scaleVal,testRotation)

    %inputs:
    %image -> image to be processed
    %backgroundImage -> background image to compare "image" to (must be
    %                                                   same size)
    %alignInfo -> struct containing parameters for alignment
    %angleGuess -> guess for rotation angle (eliminates 180 degree
    %                                                   degeneracy)
    %testRotation -> binary, true if conifirmation of the rotation angle is
    %                                                   desired

    s = size(image);
    
    if scaleVal ~= 1
        image1 = rescaleImage(image,scaleVal);
    else
        image1 = image;
    end
    
    mask = image1 >= alignInfo.bodyThreshold;
    count = 1;
    while sum(mask(:)) == 0 && count < alignInfo.bodyThreshold
        mask = image1 >= alignInfo.bodyThreshold - count;
        count = count + 1;
    end
        
    props = regionprops(mask,'Centroid','Area','PixelIdxList');
    [~,sortIdx] = sort([props.Area],'descend');
    image1_2 = uint8(zeros(size(image1)));
    if ~isempty(sortIdx)
        image1_2(props(sortIdx(1)).PixelIdxList) = image1(props(sortIdx(1)).PixelIdxList);
    end
    
    %initial alignment
    [rotationAngle,X,Y,~,~,loopImage] = ...
        alignTwoImages(backgroundImage,image1_2,angleGuess,alignInfo,image1);
    
    
    
    se = strel('square',alignInfo.erodeSize);
    midIdx = (-alignInfo.yRange:alignInfo.yRange) + round(s(1)/2);
    b = loopImage;
    c = loopImage;
    b(~alignInfo.asymImage) = 0;
    c(~fliplr(alignInfo.asymImage)) = 0;
    b = b(midIdx,:);
    sX = length(b(1,:));
    c = c(midIdx,:);
    
    b(b > alignInfo.bodyThreshold) = 0; b(b < alignInfo.minImageValue)=0;
    b = imerode(b,se);
    c(c > alignInfo.bodyThreshold) = 0; c(c < alignInfo.minImageValue)=0;
    c = imerode(c,se);
    
    sumB = sum(b(:)>0);
    sumC = sum(c(:)>0);
    
    asymValue = sumB - sumC;
    
    
    %check to see if direction is correct
    if testRotation
        
        if asymValue > 0
            
            g = mod(rotationAngle-180,360);
            [tempAngle,tempX,tempY,~,~,tempLoopImage] = ...
                alignTwoImages(backgroundImage,image1_2,g,alignInfo,image1);
            
            
            b = tempLoopImage;
            c = tempLoopImage;
            b(~alignInfo.asymImage) = 0;
            c(~fliplr(alignInfo.asymImage)) = 0;
            b = b(midIdx,:);
            sX = length(b(1,:));
            c = c(midIdx,:);
            
            b(b > alignInfo.bodyThreshold) = 0; b(b < alignInfo.minImageValue)=0;
            b = imerode(b,se);
            c(c > alignInfo.bodyThreshold) = 0; c(c < alignInfo.minImageValue)=0;
            c = imerode(c,se);

            sumB = sum(b(:)>0);
            sumC = sum(c(:)>0);

            asymValue2 = sumB - sumC;
            
            
            if asymValue2 < asymValue
                rotationAngle = tempAngle;
                X = tempX;
                Y = tempY;
                loopImage = tempLoopImage;
                asymValue = asymValue2;
            end
            
        end
        
    end