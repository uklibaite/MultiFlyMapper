function splitAvi_twoFiles(filePath,savePath,duration,ffmpegPath)

    if nargin < 4 || isempty(ffmpegPath)
        fprintf('No ffmpeg path provided. Using /opt/local/bin/\n');
        ffmpegPath = '/opt/local/bin/';
    end
    
    if nargin < 2 || isempty(savePath)
        [dir,~,~] = fileparts(filePath);
        savePath = [dir '/'];
    end
    
    if nargin < 3 || isempty(duration)
        vidObj = VideoReader(filePath);
        duration  = vidObj.Duration;
    end
    
    
    endPoint1 = num2str(ceil(duration/2));
    startPoint2 = endPoint1;
    