function [duration,status] = findMovieDuration(file,mplayerPath)

    if nargin < 2 || isempty(mplayerPath)
        mplayerPath = '/opt/local/bin/';
    end
    
    
   
    [status,duration] = unix([mplayerPath 'mplayer -identify -frames 0 -vo null -nosound ' ...
        file ' 2>&1 | awk -F= ''/LENGTH/{print $2}''']);
    
    
    duration = str2double(duration);