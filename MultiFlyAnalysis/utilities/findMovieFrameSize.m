function [imageSize,statuses] = findMovieFrameSize(filePath,tcPath)

    if nargin < 2 || isempty(tcPath)
        tcPath = '/opt/local/bin/';
    end
    
    sizes = cell(2,1);
    statuses = cell(2,1);
    
    [statuses{1},sizes{1}] = unix([tcPath 'tcprobe -i ' filePath ...
        ' 2>&1 | grep "import frame" | cut -f 2 -d "g" | cut -f 2 -d " " | cut -f 1 -d "x" ']);
    
    [statuses{2},sizes{2}] = unix([tcPath 'tcprobe -i ' filePath ...
        ' 2>&1 | grep "import frame" | cut -f 2 -d "g" | cut -f 2 -d " " | cut -f 2 -d "x" ']);
   
     
    imageSize = zeros(2,1);
    imageSize(1) = str2double(sizes{1});
    imageSize(2) = str2double(sizes{2});