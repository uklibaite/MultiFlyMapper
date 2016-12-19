function [flyData,p] = runPreprocessingEG2(filePath,p)

%create movie at 1fps, split into 300 frame segments, load segments in
    %parallel, save smaller movies and datafiles, delete movies
%create movie at 1fps, split into 300 frame segments 
% unix([p.ffmpegPath 'ffmpeg -i ' filePath ' -filter:v "crop=1220:1080:400:1" ' p.savePath 'adjusted1.avi'])
% filePath2 = [p.savePath 'adjusted1.avi'];
    
vidObj = VideoReader(filePath);
testImage = read(vidObj,1);
s = size(testImage);
p.s = s;
N = vidObj.NumberOfFrames;
frameRate = vidObj.FrameRate;

unix(['mkdir ' p.savePath '/splitMovies/']);

numSplits = floor(N/300);
splitFrames = zeros(numSplits,2);
for i = 1:numSplits
    splitFrames(i,1) = 1+(i-1)*300;
    splitFrames(i,2) = i*300;
end


boxX = -49:50;
boxY = -49:50;
boxB = -99:100;
dataCell = cell(numSplits,1);

movies = cell(1,numSplits);
for i = 1:numSplits
    movies{i} = [p.savePath sprintf('splitMovies/split_%04d.avi',i)];
end

tic
for i = 1:numSplits
    i
    images = read(vidObj,[splitFrames(i,1) splitFrames(i,2)]);
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
                
                RB1 = round(c1+50);
                LT1 = round(c1-50);
                RB2 = round(c2+50);
                LT2 = round(c2-50);
                
                [xs,ys] = boxIn(boxX+round(c1(1)),boxY+round(c1(2)),s);
                outImages(1:100,1:100,j) = images(ys,xs,j);
                
                [xs,ys] = boxIn(boxX+round(c2(1)),boxY+round(c2(2)),s);
                outImages(1:100,101:200,j) = images(ys,xs,j);
                
                dataC(j,:) = [1 2 LT1 RB1 LT2 RB2];
                
            elseif d < 200
                
                RB1 = round((c1+c2)/2+100);
                LT1 = round((c1+c2)/2-100);
                
                [xs, ys] = boxIn(boxB+round((c1(1)+c2(1))/2),boxB+round((c1(2)+c2(2))/2),s);
                outImages(:,:,j) = images(ys,xs,j);
                
                dataC(j,:) = [1 2 LT1 RB1 0 0 0 0];
                
            end
            
        elseif numBig == 1
            
            c1 = props(sortIdx(1)).Centroid;
            
            RB1 = round(c1+50);
            LT1 = round(c1-50);
            
            [xs,ys] = boxIn(boxB+round(c1(1)),boxB+round(c1(2)),s);
            outImages(:,:,j) = images(ys,xs,j);
            dataC(j,:) = [1 1 LT1 RB1 0 0 0 0];
            
        end
        
        
    end
    
    flyMovie = VideoWriter([movies{i}(1:end-4) 'p.avi']);
    open(flyMovie);
    for t = 1:size(images,3)
        writeVideo(flyMovie,uint8(outImages(:,:,t)));
    end
    close(flyMovie);
    dataCell{i} = dataC;
    
end
toc

DC = combineCells(dataCell);
[status,~] = unix([p.mencoderPath 'mencoder -oac copy -ovc copy -o ' p.savePath 'prig.avi ' p.savePath 'splitMovies/' '*p.avi']);
save([p.savePath 'prigData.mat'],'DC');


%run fly segmentation code + alignment
p.filePath = [p.savePath 'prig.avi'];
p.dataPath = [p.savePath 'prigData.mat'];
p.dataNum = 1;

% **delete split movies and adjusted movie
[status,~] = unix(['rm -r ' p.savePath 'splitMovies/']);


flyData = [];

end