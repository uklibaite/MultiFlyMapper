function  splitAvi_combine(fileName,savePath,splitLength,ffmpegPath,vidObj)

    readout = 50;

    if nargin < 3 || isempty(splitLength)
        splitLength = 300;
    else
        splitLength = ceil(splitLength);
    end
    
    if nargin < 4 || isempty(ffmpegPath)
        fprintf('No ffmpeg path provided. Using /opt/local/bin/\n');
        ffmpegPath = '/opt/local/bin/';
    end
    
    if isempty(savePath)
        [dir,name,~] = fileparts(fileName);
        dir = [dir '/' name '/splitMovies/'];
        fprintf(['No savePath, making dir ' dir '\n']);
    else
        dir = [savePath 'splitMovies/'];
        fprintf(['Saving to ' dir '\n']);
    end
        
    [status,~]=unix(['ls ' dir]);
    if status ~= 0
        [~,~] = unix(['mkdir ' dir]);
        fprintf(['New directory made: ' dir '\n']);
    else
        fprintf(['New directory not made' '\n']);
    end
    
    if nargin < 5 || isempty(vidObj)
        vidObj = VideoReader(fileName);
    end
    
    fprintf(1,'Running Splitting Script\n');
    [~,~] = unix(['python ffmpeg-split.py -f ' fileName ' -s ' num2str(splitLength) ' -p ' dir ' -m ' ffmpegPath]);
    
    
%     fprintf(1,'Finding Excess Frame Files\n');
%     
%     movies = findImagesInFolder(dir,'avi');
%     numMovies = length(movies);
%     numLeftOut = zeros(numMovies,1);
%     numFrames = zeros(numMovies,1);
%     endFrames = zeros(numMovies,1);
%     numZeros = ceil(log(numMovies)./log(10)+1e-10);
%     totalFrames = vidObj.NumberOfFrames;
%     
%     parfor i=1:numMovies
%         
%         temp = VideoReader(movies{i});
%         numFrames(i) = temp.NumberOfFrames;
%         if i < numMovies
%             numLeftOut(i) = splitLength - numFrames(i);
%         else
%             expectedNumFrames = totalFrames - (i-1)*splitLength;
%             numLeftOut(i) = expectedNumFrames - numFrames(i);
%         end
%         endFrames(i) = (i-1)*splitLength + numLeftOut(i);
%         
%     end
%     
%     fprintf(1,'Writing Excess Frame Files\n');
%     for i=1:numMovies
%         
%         if mod(i,readout) == 0
%             fprintf(1,'\t File # %6i out of %6i\n',i,numMovies);
%         end
%         
%         if numLeftOut(i) > 0
%             a = read(vidObj,[((i-1)*splitLength+1) endFrames(i)]);
%             q = num2str(i);
%             r = [repmat('0',1,numZeros-length(q)) q];
%             b = VideoWriter([dir '/resplit_' r '.avi'],'Uncompressed AVI');
%             open(b);
%             writeVideo(b,a);
%             close(b);
%         end
%         
%         clear a b q 
%         
%     end
%     
%     
%     