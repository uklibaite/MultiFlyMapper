function [fps,status] = findMovieFrameRate(file,ffmpegPath)

    if nargin < 2 || isempty(ffmpegPath)
        ffmpegPath = '/opt/local/bin/';
    end
    
    
    [status,fps] = unix([ffmpegPath 'ffmpeg -i ' file ...
        ' 2>&1 | sed -n "s/.*, \(.*\) fp.*/\1/p"']);
        
    
    fps = str2double(fps);
    