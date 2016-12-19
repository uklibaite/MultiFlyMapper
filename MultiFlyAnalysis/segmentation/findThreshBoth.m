function [threshold1 threshold2] = findThreshBoth(inputImage,parameters)
    % find body threshold using GMM
    props = regionprops(inputImage > parameters.minImageValue,'Centroid','Area','PixelIdxList');
    [currentAreas,sortIdx] = sort([props.Area],'descend');
    numBig = sum(currentAreas > parameters.minSize);
    s = size(inputImage);
    if numBig >= 2
        
        II = inputImage(inputImage>0);
        data = double(II);
        obj = gmixPlot(data,3,[],[],true,[],[],[],10);
        [~,sortIdx] = sort(obj.mu,'descend');
        
        minVal = min(data(:));
        maxVal = max(data(:));
        xx = linspace(minVal,maxVal,10000)';
        posts = posterior(obj,xx);
        
        f1 = fit(xx,posts(:,sortIdx(1))-posts(:,sortIdx(2)),'linearinterp');
        f2 = fit(xx,posts(:,sortIdx(2))-posts(:,sortIdx(3)),'linearinterp');
        
        threshold1 = fzero(f1,.5*(obj.mu(sortIdx(2))+obj.mu(sortIdx(3))));
        threshold2 = fzero(f2,.5*(obj.mu(sortIdx(3)) + obj.mu(sortIdx(2))));
        
    else
        threshold1 = 0; threshold2 = 0;
    end