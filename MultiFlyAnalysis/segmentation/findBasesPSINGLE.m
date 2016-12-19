function [basis1] = findBasesPSINGLE(inputImages,parameters)
% returns body basis images given a set of test images including 2 flies

    se1 = strel('disk',10,0);
    imageLength = parameters.imageLength;

    inputImagesR = reshape(inputImages,imageLength,imageLength,[]);
    [s1 s2 sN] = size(inputImagesR);
    imageLength = s1;
    props = cell(1,sN);
    propsAll = cell(1,sN);
    bodyAll = zeros(size(inputImagesR));
    imAll = zeros(size(inputImagesR));
    areas = zeros(1,sN);
    for i = 1:sN
        CI = inputImagesR(:,:,i);
        propsA = regionprops(CI>parameters.bodyThreshold,'Centroid','Area','PixelIdxList');
        propsAllA = regionprops(CI>parameters.minImageValue,'Centroid','Area','PixelIdxList');
        [currentAreasAll,sortIdxAll] = sort([propsAllA.Area],'descend');
        [currentAreas,sortIdx] = sort([propsA.Area],'descend');
        pBody = zeros(imageLength,imageLength); pAll = zeros(imageLength,imageLength);
        pBody(propsA(sortIdx(1)).PixelIdxList) = CI(propsA(sortIdx(1)).PixelIdxList);
        pAll(propsAllA(sortIdxAll(1)).PixelIdxList) = CI(propsAllA(sortIdxAll(1)).PixelIdxList);
        
        bodyAll(:,:,i) = pBody;
        imAll(:,:,i) = pAll;
        areas(i) = propsA(sortIdx(1)).Area;
        props{i} = propsA;
        propsAll{i} = propsAllA;
    end
    
    fprintf(1,'Generating basis images');
    % sort by areas, align 
    
    [areasSorted areasIdx] = sort(areas,'descend');
    
    b1 = parameters.backgroundImages{1};
    p = regionprops(b1>1,'Centroid');
    pC = p.Centroid;
    bgImage1 = b1(pC(2)-74:pC(2)+75,pC(1)-74:pC(1)+75);
    bgImage1(bgImage1>1)=200;
    
    
    parameters.asymImage = [zeros(150,100) ones(150,50)];
    
    alignedAll = zeros(size(inputImagesR));
    angleGuess = 0;
    for i = 1:sN
        [X,Y,rotationAngle,loopImage1,asymVal] = alignImage(inputImagesR(:,:,i),bgImage1,parameters,angleGuess,1,1);
        alignedAll(:,:,i) = loopImage1;
    end
    
    fly1Ims = alignedAll;
    fly1Med = median(fly1Ims,3);
    fly1Body = fly1Med>parameters.bodyThreshold;
    fly1NonBody = fly1Med>parameters.minImageValue & fly1Med<parameters.bodyThreshold;
    basis1 = cell(1,2);
    f1 = find(sum(fly1Body)~=0,1,'last')-15;
    cut1 = zeros(size(fly1Body)); cut1(:,1:f1) = ones(size(fly1Body,1),f1);
    cut1B = imclose(fly1Body+fly1NonBody>=1,ones(15,15));
    basis1{1} = double(fly1Body); basis1{2} = double(fly1NonBody); basis1{3} = double(cut1B);
    
    
end
