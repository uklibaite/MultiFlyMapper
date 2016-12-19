function numSplits = splitAvi_multiple(fileName,savePath,numSplits,ffmpegPath,vidObj)

    if nargin < 4 || isempty(ffmpegPath)
        fprintf('No ffmpeg path provided. Using /opt/local/bin/\n');
        ffmpegPath = '/opt/local/bin/';
    end

    if nargin < 3 || isempty(numSplits)
        numSplits = 2;
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
    
    t = vidObj.Duration;
    %nFrames = vidObj.NumberOfFrames;
    %frameRate = vidObj.FrameRate;
    
    T = ceil(t/numSplits);
    numZeros = ceil(log(numSplits)./log(10)+1e-10);
    
    readout = 50;
    
    parfor i=1:numSplits
        
        if mod(i,readout) == 0
            fprintf(1,['\t Splitting into File #' num2str(i) ' out of ' num2str(numSplits) '\n']);
        end
        
        
        q = num2str(i);
        q = [repmat('0',1,numZeros-length(q)) q];
        
        startTime = num2str((i-1)*T);
        
        
        if i < numSplits
            runLength = num2str(T);
            %endTime = num2str(i*T);
              
            %[~,~] = unix([ffmpegPath 'mencoder -ss ' startTime ' -endpos ' endTime  ...
            %    ' ' fileName ' -ovc copy -oac copy -noskip -o ' dir 'split_' q '.avi']);
            
            %fprintf(1,[ffmpegPath 'mencoder -ss ' startTime ' -endpos ' endTime  ...
            %    ' ' fileName ' -ovc copy -oac copy -noskip -o ' dir 'split_' q '.avi\n']);
            
             [~,~] = unix([ffmpegPath 'ffmpeg -r 1 -ss ' startTime ' -t ' runLength ' -i ' fileName...
                ' -vcodec copy -acodec copy -r 1 ' dir 'split_' q '.avi']);

             %[~,~] = unix([ffmpegPath 'ffmpeg -r 1 -i ' fileName '-ss ' startTime ' -t ' runLength ...
             %    ' -vcodec copy -acodec copy -r 1 ' dir 'split_' q '.avi']);
             
        else
            %endTime = num2str(t);
                        
            [~,~] = unix([ffmpegPath 'ffmpeg -r 1 -ss ' startTime ' -i ' fileName...
                ' -vcodec copy -acodec copy -r 1 ' dir 'split_' q '.avi']);
             
            %[~,~] = unix([ffmpegPath 'mencoder -ss ' startTime ' -endpos ' endTime  ...
            %    ' ' fileName ' -ovc copy -oac copy -noskip -o ' dir 'split_' q '.avi']);
                        
            %fprintf(1,[ffmpegPath 'mencoder -ss ' startTime ' -endpos ' endTime  ...
            %    ' ' fileName ' -ovc copy -oac copy -noskip -o ' dir 'split_' q '.avi\n']);
        end
        
        
    end