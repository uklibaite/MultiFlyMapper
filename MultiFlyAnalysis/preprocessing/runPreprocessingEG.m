function [flyData,p] = runPreprocessingEG(filePath,p)

%create movie at 1fps, split into 300 frame segments 
unix([p.ffmpegPath 'ffmpeg -i ' filePath ' -filter:v "crop=1220:1080:400:1" ' p.savePath 'adjusted1.avi'])
filePath2 = [p.savePath 'adjusted1.avi'];

if mod(p.splitLength,10) ~= 0
    p.splitLength = ceil(p.splitLength/10)*10;
end


split = false;
movies = findImagesInFolder(p.savePath,'mat');
if isempty(movies)
    split = true;
end


if split
    fprintf(1,'Initializing\n');
    
    fprintf(1,'\t Opening movie file\n');
    fprintf(1,['\t' filePath1 '\n']);
    vidObj = VideoReader(filePath2);
    testImage = read(vidObj,1);
    s = size(testImage);
    p.s = s;
    N = vidObj.NumberOfFrames;
    frameRate = vidObj.FrameRate;
    
    fprintf(1,'\t Adjusting move frame rate\n');
    unix([p.ffmpegPath 'ffmpeg -i ' filePath2 ' -vcodec copy -r 1 -y ' p.savePath 'adjustedMovie.avi']);
    fprintf(1,[p.ffmpegPath 'ffmpeg -i ' filePath2 ' -vcodec copy -r 1 -y ' p.savePath 'adjustedMovie.avi']);
    filePath = [p.savePath 'adjustedMovie.avi'];
    vidObj = VideoReader(filePath);
    N = vidObj.NumberOfFrames;
    
    frames = 1:N;
    
    
    frameRate = vidObj.FrameRate;
    fprintf(1,['\t Read adjusted movie, frameRate = ' frameRate '\n']);
    adjustedFrameRate = true;
    
    splitSeconds = floor(p.splitLength / frameRate);
    while mod(splitSeconds*frameRate,1) ~= 0
        splitSeconds = splitSeconds + 1;
    end
    
    fprintf(1,'Splitting File and initializing movies\n');
    unix(['mkdir ' p.savePath '/splitMovies/']);
    fprintf([p.ffmpegPath 'ffmpeg -i ' filePath ' -acodec copy -f segment -segment_time 300 -map 0 -vcodec copy -reset_timestamps 1 ' p.savePath 'splitMovies/split_%04d.avi '])
    unix([p.ffmpegPath 'ffmpeg -i ' filePath ' -acodec copy -f segment -segment_time 300 -map 0 -vcodec copy -reset_timestamps 1 ' p.savePath 'splitMovies/split_%04d.avi '])
    movies = findImagesInFolder([p.savePath '/splitMovies/'],'avi');
    numMovies = length(movies);
    
    fprintf(1,['Split into ' num2str(numMovies) ' movies']);
    
    boxX = -49:50;
    boxY = -49:50;
    boxB = -99:100;
    dataCell = cell(numMovies,1);
    
    for i = 1:numMovies
        i
        currentMovie = VideoReader(movies{i});
        images = read(currentMovie);
        images = squeeze(images(:,:,1,:));
        outImages = zeros(200,200,size(images,3));
        dataC = zeros(size(images,3),10);
        
        for j = 1:size(images,3)
                
                props = regionprops(images(:,:,j)<100,'Centroid','Area');
                [currentAreas,sortIdx] = sort([props.Area],'descend');
                numBig = sum(currentAreas > 100);
                
                
                if numBig > 1
                    
                    c1 = props(sortIdx(1)).Centroid;
                    c2 = props(sortIdx(2)).Centroid;
                    d = sqrt(sum((c1-c2).^2));
                    
                    if d > 200
                        
                        RB1 = round(c1+75);
                        LT1 = round(c1-75);
                        RB2 = round(c2+75);
                        LT2 = round(c2-75);
                        
                        [xs,ys] = boxIn(boxX+round(c1(1)),boxY+round(c1(2)),s);
                        outImages(1:150,1:150,j) = images(ys,xs,j);
                        
                        [xs,ys] = boxIn(boxX+round(c2(1)),boxY+round(c2(2)),s);
                        outImages(1:150,151:300,j) = images(ys,xs,j);
                        
                        dataC(j,:) = [1 2 LT1 RB1 LT2 RB2];
                        
                    elseif d < 200
                        
                        RB1 = round((c1+c2)/2+150);
                        LT1 = round((c1+c2)/2-150);
                        
                        [xs, ys] = boxIn(boxB+round((c1(1)+c2(1))/2),boxB+round((c1(2)+c2(2))/2),s);
                        outImages(:,:,j) = images(ys,xs,j);
                        
                        dataC(j,:) = [1 2 LT1 RB1 0 0 0 0];
                        
                    end
                    
                elseif numBig == 1
                    
                    c1 = props(sortIdx(1)).Centroid;
                    
                    RB1 = round(c1+75);
                    LT1 = round(c1-75);
                    
                    [xs,ys] = boxIn(boxB+round(c1(1)),boxB+round(c1(2)),s);
                    outImages(:,:,j) = images(ys,xs,j);
                    dataC(j,:) = [1 1 LT1 RB1 0 0 0 0];
                    
                end
                
                
            end
        
        
    end
        
        
    
else   
end





%open save paths

% parameters adjustments
p.minImageValue = 120;
minImageValue = p.minImageValue;
p.minSize = 300;
minSize = p.minSize;
imageSize = 100;
imageLength = 100;
p.imageLength = imageLength;
L = length(firstMovie:lastMovie);
flyData1_centroid = cell(L,1);
flyData1_area = cell(L,1);
flyData1_leftUpperCorner = cell(L,1);
flyData2_centroid = cell(L,1);
flyData2_area = cell(L,1);
flyData2_leftUpperCorner = cell(L,1);
numInImages = cell(L,1);
%imageLength = 2*imageSize + 1;


% set threshold for all images for coarse tracking
bodyThreshold = 200;
backgroundT = p.minImageValue;






end