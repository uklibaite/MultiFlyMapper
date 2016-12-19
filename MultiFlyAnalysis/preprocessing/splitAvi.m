function N = splitAvi(fileName,savePath,splitLength)

    %split length is in minutes
    splitLength = round(splitLength);
    
    %[dir,name,~] = fileparts(fileName);
    %dir = [dir '/' name '/splitMovies/'];
    dir = [savePath '/splitMovies/'];
    
    [status,~]=unix(['ls ' dir]);
    if status ~= 0
        unix(['mkdir ' dir]);
    end
    
    vidObj = VideoReader(fileName);
    t = vidObj.Duration / 60;
    N = ceil(t / splitLength);
    numZeros = ceil(log(N)./log(10)+1e-10);
    
    for i=1:N
        
        if mod(i,50) == 0
            fprintf(1,'\t Splitting Movie #%5i out of %5i\n',i,N);
        end
        
        q = num2str(i);
        q = [repmat('0',1,numZeros-length(q)) q];
        
        sL = num2str(splitLength);
        sL = [repmat('0',1,2-length(sL)) sL];
        
        startTime = (i-1)*splitLength;
        hours = num2str(floor(startTime / 60));
        hours = [repmat('0',1,2-length(hours)) hours];
        minutes = num2str(mod(startTime,60));
        minutes = [repmat('0',1,2-length(minutes)) minutes];
        
        
        fprintf(1,[hours ':' minutes ':00.0 -t 00:' sL ':0.0\n']);
        
        if i < N
  
            [~,~] = unix(['/opt/local/bin/ffmpeg -r 1 -ss ' hours ':' minutes ':00.0 -t 00:' sL ':0.0 -i ' fileName...
                ' -vcodec copy -acodec copy -r 1 ' dir 'split_' q '.avi']);
            
        else
      
            [~,~] = unix(['/opt/local/bin/ffmpeg -r 1 -ss  ' hours ':' minutes ':00.0 -i ' fileName...
                ' -vcodec copy -acodec copy -r 1 ' dir 'split_' q '.avi']);
            
        end
        
    end
