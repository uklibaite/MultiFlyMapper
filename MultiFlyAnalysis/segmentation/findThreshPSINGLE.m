function [threshold] = findThreshPSINGLE(inputImage,parameters)
    % find body threshold using GMM
    props = regionprops(inputImage > parameters.minImageValue,'Centroid','Area','PixelIdxList');
    [currentAreas,sortIdx] = sort([props.Area],'descend');
    numBig = sum(currentAreas > parameters.minSize);
    s = size(inputImage);
    
    ix = s(2);
    iy = s(1);
    
    cy = props(sortIdx(1)).Centroid(2);
    cx = props(sortIdx(1)).Centroid(1);
    
    
    r = 60;
    [x, y] = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    m = ((x.^2+y.^2)<=r^2);
    
    
    cmask = zeros(iy,ix);
    cmask(m>=1) = 1;
    mIm = uint8(cmask).*inputImage;
    
    II = mIm(mIm>0);
    threshold = autoFindThreshold_gmm(II,3);
    
