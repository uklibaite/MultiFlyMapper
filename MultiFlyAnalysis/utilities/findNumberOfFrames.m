function [N,status] = findNumberOfFrames(file,tcPath)

    if nargin < 2 || isempty(tcPath)
        tcPath = '/opt/local/bin/';
    end
    
    
    [status,N] = unix([tcPath 'tcprobe -i ' file ...
        ' 2>&1 | grep "length:" | cut -f 1 -d "," | cut -f 2 -d ":" | cut -f 2 -d " " ']);
   
        
    N = str2double(N);