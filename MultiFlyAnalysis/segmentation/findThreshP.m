function [threshold] = findThreshP(inputImage,parameters)
    % find body threshold using GMM
    props = regionprops(inputImage > parameters.minImageValue,'Centroid','Area','PixelIdxList');
    [currentAreas,sortIdx] = sort([props.Area],'descend');
    numBig = sum(currentAreas > parameters.minSize);
    s = size(inputImage);
    if numBig >= 2
        ix = s(2);
        iy = s(1);
        
        cy = props(sortIdx(1)).Centroid(2);
        cx = props(sortIdx(1)).Centroid(1);
        cy2 = props(sortIdx(2)).Centroid(2);
        cx2 = props(sortIdx(2)).Centroid(1);
        
        r = 60;
        [x, y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
        m = ((x.^2+y.^2)<=r^2);
        
        [x2, y2] = meshgrid(-(cx2-1):(ix-cx2),-(cy2-1):(iy-cy2));
        m2 = ((x2.^2+y2.^2)<=r^2);
        
        cmask = zeros(iy,ix);
        cmask(m+m2>=1) = 1;
        mIm = uint8(cmask).*inputImage;
        
        II = mIm(mIm>0);
        threshold = autoFindThreshold_gmm(II,3);
    else
        threshold = 0;
    end

