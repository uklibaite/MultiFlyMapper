function [flyData,runData,p] = multiFlySegmentAndAlignSINGLE(startFrame,endFrame,p)

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
bigMovie = read(vidObj,[startFrame endFrame]);
toc

bigMovie = squeeze(bigMovie(:,:,1,:));

% read data for movie, figure out connectivities
load(p.dataPath);
CamData = DC;

% find threshold for all images
xB = randsample((1:N),300);
testImages = uint8(zeros(p.imageLength,p.imageLength,length(xB)));
for t = 1:length(xB)
    testImages(:,:,t) = bigMovie(1:p.imageLength,:,xB(t));
end


T = zeros(length(xB),1);
for t = 1:length(xB)
    if mod(t,20)==0
        t
    end
    currentImage = imcomplement(testImages(:,:,t));
    T(t) = findThreshPSINGLE(currentImage,p);
end

bodyThreshold = median(T(T~=0));
if isnan(bodyThreshold)
    bodyThreshold = p.defaultBodyThreshold;
end
fprintf('body Threshold = %i\n',bodyThreshold);
p.bodyThreshold = bodyThreshold;
 

p.asymImage = [zeros(150,100) ones(150,50)];
[basis1] = findBasesPSINGLE(imcomplement(testImages),p);

backgroundT = p.minImageValue;
%runData.frames = frames;
%runData.medianImage = medianImage;
runData.bodyThreshold = bodyThreshold;
runData.filePath = p.filePath;
runData.savePath = p.savePath;
%runData.segmentationFrames = frames;
clear numBig

centroids1 = zeros(N,2);
areas1 = zeros(N,1);
thetas1 = zeros(N,1);
shift1 = zeros(N,2);
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
    dataLine = CamData(i,:);
    cent1 = [(dataLine(3)+dataLine(5))/2 (dataLine(4)+dataLine(6))/2];
    outImage = zeros(p.imageLength,p.imageLength);
    if i == 1
        initialPhi1 = 0;
    else initialPhi1 = thetas1(i-1);
    end
    
    p1All = regionprops(CI(1:p.imageLength,1:p.imageLength)>backgroundT,'Area','PixelIdxList');
    [pA p1Idx] = sort([p1All.Area],'descend');
    
    pixIm1 = zeros(p.imageLength,p.imageLength);
    C1 = CI(1:p.imageLength,1:p.imageLength);
    pixIm1(p1All(p1Idx(1)).PixelIdxList)= C1(p1All(p1Idx(1)).PixelIdxList);
    areas1(i) = sum(sum(pixIm1>bodyThreshold));
    
    
    [X1,Y1,rotationAngle1,loopImage1,asymVal1] = alignImage(pixIm1,basis1{1},p,initialPhi1,1,testRotation);
    shift1(i,:) = [X1,Y1];
    centroids1(i,:) = cent1+[X1,Y1];
    outImage(1:p.imageLength,1:p.imageLength) = uint8(loopImage1);
    thetas1(i) = rotationAngle1;
    NII(i) = 0;
    outImages(:,:,i) = outImage;
    
    if mod(i,1000)==0
       fprintf(1,['Tracking will continue == ' num2str(continueTrack) '\n']);
    end
    
end


p.frames = 1:N;
basis1Area = sum(sum(basis1{1}));
sdiff1 = p.rescaleBaseline-basis1Area; sig1 = sign(sdiff1); sfrac1 = sqrt(abs(sdiff1/p.rescaleBaseline)*100);
scaleVal1 = 1+(sig1*sfrac1*2)/100;
runData.scaleVal1 = scaleVal1;


flyMovie1 = VideoWriter([p.savePath 'flies_aligned_SINGLE' ]);
open(flyMovie1);
for i = p.frames
    if mod(i,1000)==0
        i
    end
    writeVideo(flyMovie1,uint8(rescaleImage(outImages(1:p.imageLength,1:p.imageLength,i),scaleVal1)));
end
close(flyMovie1);

flyData{1}.centroid = centroids1(p.frames,:);
flyData{1}.area = areas1(p.frames,:);
flyData{1}.shift = shift1(p.frames,:);
flyData{1}.theta = thetas1(p.frames,:);

runData.fly1Area = median(flyData{1}.area);
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
