function [flyData,runData,p] = multiFlySegmentAndAlign(startFrame,endFrame,p)

fprintf(1,'Initializing Alignment\n');

fprintf(1,'\t Opening movie file\n');
vidObj = VideoReader(p.filePath);
testImage = read(vidObj,1);
s = size(testImage);
p.s = s;

ix = s(2);
iy = s(1);


N = vidObj.NumberOfFrames;
if endFrame <= N
    N = endFrame;
end

frameRate = vidObj.FrameRate;

tic
bigMovie = read(vidObj,[startFrame N]);
toc

bigMovie = squeeze(bigMovie(:,:,1,:));

% read data for movie, figure out connectivities
load(p.dataPath);
CamData = DC;

CamZeros = sum(CamData == 0,2);
CamIsBoth0 = CamZeros == 0;
CamIsBoth4 = CamZeros>=4;
CamIsBoth8 = CamZeros>=8;
Separate = find(CamIsBoth0==1);

% find threshold for all images
xB = randsample(Separate(Separate<N),300);
testImages = uint8(zeros(p.imageLength,p.imageLength*2,length(xB)));
for t = 1:length(xB)
    testImages(:,:,t) = bigMovie(1:p.imageLength,:,xB(t));
end


T = zeros(length(xB),1);
for t = 1:length(xB)
    if mod(t,20)==0
        t
    end
    currentImage = imcomplement(testImages(:,:,t));
    T(t) = findThreshP(currentImage,p);
end

bodyThreshold = median(T(T~=0));
if isnan(bodyThreshold)
    bodyThreshold = p.defaultBodyThreshold;
end
fprintf('body Threshold = %i\n',bodyThreshold);
p.bodyThreshold = bodyThreshold;
 

p.asymImage = [zeros(150,100) ones(150,50)];
[basis1 basis2] = findBasesP(imcomplement(testImages),p);

backgroundT = p.minImageValue;
%runData.frames = frames;
%runData.medianImage = medianImage;
runData.bodyThreshold = bodyThreshold;
runData.filePath = p.filePath;
runData.savePath = p.savePath;
%runData.segmentationFrames = frames;
clear numBig

centroids1 = zeros(N,2);
centroids2 = zeros(N,2);

areas1 = zeros(N,1);
areas2 = zeros(N,1);

thetas1 = zeros(N,1);
thetas2 = zeros(N,1);

shift1 = zeros(N,2);
shift2 = zeros(N,2);

NII = zeros(N,1);
FP = zeros(N,1);
testRotation = 1;
continueTrack = 1;
copThresh = 10;
i = 0;
Ncap = N-mod(N,1000);



while continueTrack ~= 0 && i <= Ncap-1
    i = i+1;
    if mod(i,1000)==0
        i
    end
    CI = imcomplement(bigMovie(:,:,i));
    isBoth = CamIsBoth0(i,1);
    dataLine = CamData(i,:);
    cent1 = [(dataLine(3)+dataLine(5))/2 (dataLine(4)+dataLine(6))/2];
    cent2 = [(dataLine(7)+dataLine(9))/2 (dataLine(8)+dataLine(10))/2];
    outImage = zeros(p.imageLength,2*p.imageLength);
    if i == 1
        initialPhi1 = 0; initialPhi2 = 0;
    else initialPhi1 = thetas1(i-1); initialPhi2 = thetas2(i-1);
    end
    
    if isBoth
        
        p1All = regionprops(CI(1:p.imageLength,1:p.imageLength)>backgroundT,'Area','PixelIdxList');
        p2All = regionprops(CI(1:p.imageLength,p.imageLength+(1:p.imageLength))>backgroundT,'Area','PixelIdxList');
        [pA p1Idx] = sort([p1All.Area],'descend'); [pA p2Idx] = sort([p2All.Area],'descend');
        
        pixIm1 = zeros(p.imageLength,p.imageLength); pixIm2 = zeros(p.imageLength,p.imageLength);
        C1 = CI(1:p.imageLength,1:p.imageLength); C2 = CI(1:p.imageLength,p.imageLength+(1:p.imageLength));
        
        pixIm1(p1All(p1Idx(1)).PixelIdxList)= C1(p1All(p1Idx(1)).PixelIdxList);
        pixIm2(p2All(p2Idx(1)).PixelIdxList)= C2(p2All(p2Idx(1)).PixelIdxList);
        
        areas1(i) = sum(sum(pixIm1>bodyThreshold));
        areas2(i) = sum(sum(pixIm2>bodyThreshold));
        
        
        if areas1(i)<areas2(i)
            
            tempArea = areas1(i);
            areas1(i) = areas2(i);
            areas2(i) = tempArea;
            centroids1(i,:) = cent2;
            centroids2(i,:) = cent1;
            tempPixIm = pixIm1;
            pixIm1 = pixIm2;
            pixIm2 = tempPixIm;
            
        end
        
        [X1,Y1,rotationAngle1,loopImage1,asymVal1] = alignImage(pixIm1,basis1{1},p,initialPhi1,1,testRotation);
        [X2,Y2,rotationAngle2,loopImage2,asymVal2] = alignImage(pixIm2,basis2{1},p,initialPhi2,1,testRotation);
        
        shift1(i,:) = [X1,Y1];
        shift2(i,:) = [X2,Y2];
        
        centroids1(i,:) = cent1+[X1,Y1];
        centroids2(i,:) = cent2+[X1,Y1];
        
        outImage(1:p.imageLength,1:p.imageLength) = uint8(loopImage1);
        outImage(1:p.imageLength,p.imageLength+(1:p.imageLength)) = uint8(loopImage2);
        
        thetas1(i) = rotationAngle1;
        thetas2(i) = rotationAngle2;
        
        NII(i) = 2;
        
    else
        
        if cent1(1) >= 0
            ULCTemp = [dataLine(3) dataLine(4)];
        else
            ULCTemp = [dataLine(7) dataLine(8)];
        end
        props = regionprops(CI>bodyThreshold,'Centroid','Area','PixelIdxList');
        propsAll = regionprops(CI>backgroundT,'Area','PixelIdxList');
        [currentAreas,sortIdx] = sort([props.Area],'descend');
        [currentAreasAll,sortIdxAll] = sort([propsAll.Area],'descend');
        numBig = sum(currentAreasAll>p.minSize);
        numBodies = sum(currentAreas>p.minSize);
        
        if numBig >= 2 && numBodies >= 2
            % have to fix for original centroid of image (when
            % ~isBoth)
            
            centroids1(i,:) = props(sortIdx(1)).Centroid;
            centroids2(i,:) = props(sortIdx(2)).Centroid;
            areas1(i) = props(sortIdx(1)).Area;
            areas2(i) = props(sortIdx(2)).Area;
            xs = (1-p.imageSize:p.imageSize) + round(centroids1(i,1));
            ys = (1-p.imageSize:p.imageSize) + round(centroids1(i,2));
            [xs,ys] = boxIn(xs,ys,s);
            C = CI(ys,xs);
            propsTemp = regionprops(C>p.minImageValue,'Centroid','Area','PixelIdxList');
            [currentAreasTemp,sortIdxTemp] = sort([propsTemp.Area],'descend');
            
            pixIm = zeros(p.imageLength,p.imageLength);
            pixIm(propsTemp(sortIdxTemp(1)).PixelIdxList) = C(propsTemp(sortIdxTemp(1)).PixelIdxList);
            [X1,Y1,rotationAngle1,loopImage1,asymVal1] = alignImage(pixIm,basis1{1},p,initialPhi1,1,testRotation);
            outImage(1:p.imageLength,1:p.imageLength) = uint8(loopImage1);
            
            xs = (1-p.imageSize:p.imageSize) + round(centroids2(i,1));
            ys = (1-p.imageSize:p.imageSize) + round(centroids2(i,2));
            [xs,ys] = boxIn(xs,ys,s);
            C = CI(ys,xs);
            propsTemp = regionprops(CI(ys,xs)>p.minImageValue,'Centroid','Area','PixelIdxList');
               
            [currentAreasTemp,sortIdxTemp] = sort([propsTemp.Area],'descend');
            if sum(currentAreasTemp>p.minSize)>=2
                distCheck = zeros(1,length(propsTemp));
                for dd = 1:length(propsTemp)
                distCheck(dd) = sqrt(sum((propsTemp(dd).Centroid-[75 75]).^2));
                end
                tempMinDistId = find(distCheck == min(distCheck));
                if propsTemp(tempMinDistId).Area > p.minSize
                sortIdxTemp = tempMinDistId;
                end
            end
            
            
            pixIm = zeros(p.imageLength,p.imageLength);
            pixIm(propsTemp(sortIdxTemp(1)).PixelIdxList) = C(propsTemp(sortIdxTemp(1)).PixelIdxList);
            [X2,Y2,rotationAngle2,loopImage2,asymVal2] = alignImage(pixIm,basis2{1},p,initialPhi2,1,testRotation);
            outImage(1:p.imageLength,p.imageLength+(1:p.imageLength)) = uint8(loopImage2);
            NII(i) = 2;
            FP(i) = 1;
            shift1(i,:) = [X1,Y1];
            shift2(i,:) = [X2,Y2];
            centroids1(i,:) = ULCTemp+centroids1(i,:)+[X1,Y1];
            centroids2(i,:) = ULCTemp+centroids2(i,:)+[X2,Y2];
            thetas1(i) = rotationAngle1;
            thetas2(i) = rotationAngle2;
            
            
        elseif numBig == 1 && size(props,1)>1 && numBodies >= 2
            
            tempTwo = zeros(size(CI));
            tempTwoB = zeros(size(CI));
            tempTwo(propsAll(sortIdxAll(1)).PixelIdxList)=1;
            tempTwoB(props(sortIdx(1)).PixelIdxList) = 1;
            tempTwoB(props(sortIdx(2)).PixelIdxList) = 1;
            tempTwoB = imopen(tempTwoB,ones(5,5));
            
            tempTwo = zeros(size(CI));
            tempTwo(propsAll(sortIdxAll(1)).PixelIdxList)=1;
            BW = im2bw(tempTwo,.5);
            
            tempBody = zeros(size(CI));
            tempBody(props(sortIdx(1)).PixelIdxList)=1;
            tempBody(props(sortIdx(2)).PixelIdxList)=2;
            D1 = -bwdist(tempBody==1);
            D2 = -bwdist(tempBody==2);
            D1(~BW)=-Inf; D2(~BW)=-Inf;
            mask1 = D1>=D2; mask1(~BW)=0; m1props = regionprops(mask1,'Centroid','Area','PixelIdxList');
            mask2 = D2>=D1; mask2(~BW)=0; m2props = regionprops(mask2,'Centroid','Area','PixelIdxList');
            Cm1 = CI.*uint8(mask1);
            Cm2 = CI.*uint8(mask2);
            
            % set centroids and areas using masks
            centroids1(i,:) = props(sortIdx(1)).Centroid;
            centroids2(i,:) = props(sortIdx(2)).Centroid;
            areas1(i) = props(sortIdx(1)).Area;
            areas2(i) = props(sortIdx(2)).Area;
            
            % first pixIm
            xs = (1-p.imageSize:p.imageSize) + round(centroids1(i,1));
            ys = (1-p.imageSize:p.imageSize) + round(centroids1(i,2));
            [xs,ys] = boxIn(xs,ys,s);
            Cm1 = Cm1(ys,xs);
            
            % second pixIm
            xs = (1-p.imageSize:p.imageSize) + round(centroids2(i,1));
            ys = (1-p.imageSize:p.imageSize) + round(centroids2(i,2));
            [xs,ys] = boxIn(xs,ys,s);
            Cm2 = Cm2(ys,xs);
            
            p1All = regionprops(Cm1>backgroundT,'Area','PixelIdxList');
            p2All = regionprops(Cm2>backgroundT,'Area','PixelIdxList');
            [pA p1Idx] = sort([p1All.Area],'descend'); [pA p2Idx] = sort([p2All.Area],'descend');
            
            
            pixIm1 = zeros(p.imageLength,p.imageLength); pixIm2 = zeros(p.imageLength,p.imageLength);
            pixIm1(p1All(p1Idx(1)).PixelIdxList)= Cm1(p1All(p1Idx(1)).PixelIdxList);
            pixIm2(p2All(p2Idx(1)).PixelIdxList)= Cm2(p2All(p2Idx(1)).PixelIdxList);
            
            areas1(i) = sum(sum(pixIm1>bodyThreshold));
            areas2(i) = sum(sum(pixIm2>bodyThreshold));
            
            [X1,Y1,rotationAngle1,loopImage1,asymVal1] = alignImage(pixIm1,basis1{1},p,initialPhi1,1,testRotation);
            outImage(1:p.imageLength,1:p.imageLength) = uint8(loopImage1.*basis1{3});
            [X2,Y2,rotationAngle2,loopImage2,asymVal2] = alignImage(pixIm2,basis2{1},p,initialPhi2,1,testRotation);
            outImage(1:p.imageLength,p.imageLength+(1:p.imageLength)) = uint8(loopImage2.*basis2{3});
            
            shift1(i,:) = [X1,Y1];
            shift2(i,:) = [X2,Y2];
            centroids1(i,:) = ULCTemp+centroids1(i,:)+[X1 Y1];
            centroids2(i,:) = ULCTemp+centroids2(i,:)+[X2 Y2];
            thetas1(i) = rotationAngle1;
            thetas2(i) = rotationAngle2;
            
            NII(i)=2;
            
        else
            sig = .7;
            tempTwo = zeros(size(CI));
            tempTwo(propsAll(sortIdxAll(1)).PixelIdxList)=1;
            BW = im2bw(tempTwo,.5);
            D = -bwdist(~BW);
            D(~BW) = -Inf;
            D2=imhmin(D,sig);
            L = watershed(D2);
            
            while max(max(L))>=4 && sig<=1
                sig = sig+.05;
                D2 = imhmin(D,sig);
                L = watershed(D2);
            end
            
            sumX = [];
            for k = 1:max(max(L))
                sumX(k) = sum(sum(L==k));
            end
            
            [sortedSums sumIdx] = sort(sumX,'descend');
            if sum(sortedSums(3:end))>=1000
                NewL = zeros(size(L));
                NewL(L==sumIdx(2)) = 1;
                for y = 3:length(sortedSums)
                    NewL(L==sumIdx(y)) = 2;
                end
                mask1 = zeros(size(CI));
                mask1 = imclose((NewL==1),p.openingFilter);
                mask2 = zeros(size(CI)); mask2 = imclose((NewL==2),p.openingFilter);
            else
                mask1 = zeros(size(CI));
                mask1 = imclose(L==sumIdx(2),p.openingFilter);
                mask2 = zeros(size(CI));
            end
            
            
            Cm1 = CI.*uint8(mask1);
            Cm2 = CI.*uint8(mask2);
            props1 = regionprops(Cm1>=bodyThreshold,'Area','Centroid','PixelIdxList');
            props2 = regionprops(Cm2>=bodyThreshold,'Area','Centroid','PixelIdxList');
            [t idx1] = sort([props1.Area],'descend');
            [t2 idx2] = sort([props2.Area],'descend');
            
            
            
            if ~isempty(idx1)&&~isempty(idx2)
                
                centroids1(i,:) = props1(idx1).Centroid;
                centroids2(i,:) = props2(idx2).Centroid;
                
                % first pixIm
                xs = (1-p.imageSize:p.imageSize) + round(centroids1(i,1));
                ys = (1-p.imageSize:p.imageSize) + round(centroids1(i,2));
                [xs,ys] = boxIn(xs,ys,s);
                Cm1 = Cm1(ys,xs);
                
                % second pixIm
                xs = (1-p.imageSize:p.imageSize) + round(centroids2(i,1));
                ys = (1-p.imageSize:p.imageSize) + round(centroids2(i,2));
                [xs,ys] = boxIn(xs,ys,s);
                Cm2 = Cm2(ys,xs);
                
                
                p1All = regionprops(Cm1>backgroundT,'Area','PixelIdxList');
                p2All = regionprops(Cm2>backgroundT,'Area','PixelIdxList');
                [pA p1Idx] = sort([p1All.Area],'descend'); [pA p2Idx] = sort([p2All.Area],'descend');
                
                pixIm1 = zeros(p.imageLength,p.imageLength); pixIm2 = zeros(p.imageLength,p.imageLength);
                pixIm1(p1All(p1Idx(1)).PixelIdxList)= Cm1(p1All(p1Idx(1)).PixelIdxList);
                pixIm2(p2All(p2Idx(1)).PixelIdxList)= Cm2(p2All(p2Idx(1)).PixelIdxList);
                
                areas1(i) = sum(sum(pixIm1>bodyThreshold));
                areas2(i) = sum(sum(pixIm2>bodyThreshold));
                
                [X1,Y1,rotationAngle1,loopImage1,asymVal1] = alignImage(pixIm1,basis1{1},p,initialPhi1,1,testRotation);
                outImage(1:p.imageLength,1:p.imageLength) = uint8(loopImage1.*basis1{3});
                [X2,Y2,rotationAngle2,loopImage2,asymVal2] = alignImage(pixIm2,basis2{1},p,initialPhi2,1,testRotation);
                outImage(1:p.imageLength,p.imageLength+(1:p.imageLength)) = uint8(loopImage2.*basis2{3});
                
                shift1(i,:) = [X1,Y1];
                shift2(i,:) = [X2,Y2];
                centroids1(i,:) = ULCTemp+centroids1(i,:)+[X1 Y1];
                centroids2(i,:) = ULCTemp+centroids2(i,:)+[X2 Y2];
                thetas1(i) = rotationAngle1;
                thetas2(i) = rotationAngle2;
                
                NII(i) = 2;
                
            elseif ~isempty(idx1)
                centroids1(i,:) = props1(idx1(1)).Centroid;
                if i == 1
                    centroids2(i,:) = centroids1(i,:);
                else
                    centroids2(i,:) = centroids2(i-1,:);
                end
                areas1(i) = props1(idx1(1)).Area;
                areas2(i) = 0;
                xs = (1-p.imageSize:p.imageSize) + round(centroids1(i,1));
                ys = (1-p.imageSize:p.imageSize) + round(centroids1(i,2));
                [xs,ys] = boxIn(xs,ys,s);
                Ct1 = Cm1(ys,xs);
                outImage(1:p.imageLength,1:p.imageLength) = Ct1;
                
                NII(i) = 1;
                
            else
                NII(i) = 0;
                
            end
        end
    end
    outImages(:,:,i) = outImage;
    
    if mod(i,1000)==0
        
        stillCourt = medfilt1(double(NII(1:i)~=2),250);
        B = bwconncomp(stillCourt);
        C = returnCellLengths(B.PixelIdxList);
        Idx = find(C>copThresh,1,'first');
        
        if ~isempty(Idx)
            continueTrack = 0;
        else
            fprintf(1,['Tracking will continue == ' num2str(continueTrack) '\n']);
        end
        
        
    end
    
end


if ~isempty(Idx)
    frames = 1:(B.PixelIdxList{Idx}(1)-1);
else
    frames = 1:100000;
end

% p.frames = frames(1:end-500);
p.frames = frames;

basis1Area = sum(sum(basis1{1}));
basis2Area = sum(sum(basis2{1}));
sdiff1 = p.rescaleBaseline-basis1Area; sig1 = sign(sdiff1); sfrac1 = sqrt(abs(sdiff1/p.rescaleBaseline)*100);
sdiff2 = p.rescaleBaseline-basis2Area; sig2 = sign(sdiff2); sfrac2 = sqrt(abs(sdiff2/p.rescaleBaseline)*100);
scaleVal1 = 1+(sig1*sfrac1*2)/100;
scaleVal2 = 1+(sig2*sfrac2*2)/100;
runData.scaleVal1 = scaleVal1;
runData.scaleVal2 = scaleVal2;


flyMovie1 = VideoWriter([p.savePath 'flies_aligned_female' ]);
flyMovie2 = VideoWriter([p.savePath 'flies_aligned_male' ]);
open(flyMovie1);
open(flyMovie2);
for i = p.frames
    if mod(i,1000)==0
        i
    end
    writeVideo(flyMovie1,uint8(rescaleImage(outImages(1:p.imageLength,1:p.imageLength,i),scaleVal1)));
    writeVideo(flyMovie2,uint8(rescaleImage(outImages(1:p.imageLength,p.imageLength+(1:p.imageLength),i),scaleVal2)));
end
close(flyMovie1);
close(flyMovie2);


flyData{1}.centroid = centroids1(p.frames,:);
flyData{1}.area = areas1(p.frames,:);
flyData{1}.shift = shift1(p.frames,:);
flyData{1}.theta = thetas1(p.frames,:);

flyData{2}.centroid = centroids2(p.frames,:);
flyData{2}.area = areas2(p.frames,:);
flyData{2}.shift = shift2(p.frames,:);
flyData{2}.theta = thetas2(p.frames,:);
runData.fly1Area = median(flyData{1}.area(NII(p.frames)==2));
runData.fly2Area = median(flyData{2}.area(NII(p.frames)==2));
p.flyData = flyData;


runData.framesWithTwo = NII(p.frames)==2;
runData.numberInImages = NII;
end



function [xs,ys] = boxIn(xs,ys,s)

if xs(1) < 1
    xs = xs - min(xs) + 1;
else if xs(end) > s(2)
        xs = xs - (xs(end) - s(2));
    end
end

if ys(1) < 1
    ys = ys - min(ys) + 1;
else if ys(end) > s(1)
        ys = ys - (ys(end) - s(1));
    end
end

end

