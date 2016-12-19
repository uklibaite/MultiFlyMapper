function image = makeEdgesZero(image,edgeThreshold)
    
    image(1:edgeThreshold,:) = 0;
    image(:,1:edgeThreshold) = 0;
    image(end-edgeThreshold:end,:) = 0;
    image(:,end-edgeThreshold:end) = 0;

    