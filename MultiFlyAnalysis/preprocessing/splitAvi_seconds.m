function [N,splitMovies] = splitAvi_seconds(fileName,savePath,splitLength,ffmpegPath,vidObj)

    if nargin < 4 || isempty(ffmpegPath)
        fprintf('No ffmpeg path provided. Using /opt/local/bin/\n');
        ffmpegPath = '/opt/local/bin/';
    end

    %split length is in minutes
    splitLength = round(splitLength*100)/100;
    if splitLength <= 0
        splitLength = 1;
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
    
    %if isnumeric(vidObj)
    %        t = vidObj;
    %    else
    %        t = vidObj.Duration;
    %end
    %N = ceil(t / splitLength);
    %numZeros = ceil(log(N)./log(10)+1e-10);
    
    if nargin < 5 || isempty(vidObj)
        vidObj = VideoReader(fileName);
    end
    t = vidObj.Duration;
    N = ceil(t / splitLength);
    numZeros = ceil(log(N)./log(10)+1e-10);
    
    
    readout = 50;
    
    if nargout > 1
        splitMovies = cell(N,1);
        makeSplitMovies = true;
    else
        makeSplitMovies = false;
    end
    
    parfor i=1:N
        
        if mod(i,readout) == 0
            fprintf(1,['\t Splitting into File #' num2str(i) ' out of ' num2str(N) '\n']);
        end
            
        q = num2str(i);
        q = [repmat('0',1,numZeros-length(q)) q];
        
        startTime = num2str((i-1)*splitLength);
        sL = num2str(splitLength);
        
        
        if makeSplitMovies
            splitMovies{i} = [dir 'split_' q '.avi'];
        end
        
        
        if i < N
            
            [~,~] = unix([ffmpegPath 'ffmpeg -r 1 -i ' fileName ' -ss ' startTime ' -t ' sL ...
                ' -vcodec copy -r 1 -y ' dir 'split_' q '.avi']);

        else
            
             [~,~] = unix([ffmpegPath 'ffmpeg -r 1 -i ' fileName ' -ss ' startTime ...
                ' -vcodec copy -r 1 -y ' dir 'split_' q '.avi']);
            
        end
        
    end
    
