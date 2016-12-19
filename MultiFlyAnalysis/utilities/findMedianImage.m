function medianImage = findMedianImage(vidObj,N)

    N = min(N,vidObj.NumberOfFrames);

    if vidObj.NumberOfFrames <= N
        
        images = read(vidObj,[1 vidObj.NumberOfFrames]);
    
    else
    
        image = read(vidObj,randi(vidObj.NumberOfFrames));
        s = size(image);
    
        images = zeros(s(1),s(2),N);
        images(:,:,1) = image(:,:,1);
        
        for i=2:N
            image = read(vidObj,randi(vidObj.NumberOfFrames));
            images(:,:,i) = image(:,:,1);
        end
        
    end
    
    medianImage = median(images,3);