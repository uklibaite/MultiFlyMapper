function idx = makeBradyMovie_avis_readers(file_path,readers,groups,subX,subY,complementImage)

    %inputs:

    %file_path  --  where the images should be saved
    %movies     --  Lx1 cell array with strings pointing to the movies
    %groups     --  Nx3 array, first column is the movie #, second and third
    %               are the start and end frames
    %subX,subY  --  dimensions of the output movie (in # of flies)

    
    if nargin < 6 || isempty(complementImage)
        complementImage = false;
    end
    
    
    N = length(groups(:,1));
    if N > subX*subY
        N = subX*subY;
        idx = randperm(length(groups(:,1)),N);
        groups = groups(idx,:);
    else
        idx = 1:N;
        subX = ceil(sqrt(N));
        subY = subX;
        %subY = ceil(N / subX);
    end
    
    
    lengths = groups(:,3) - groups(:,2) + 1;
    
    runLength = max(lengths);
    numZeros = ceil(log(runLength)./log(10));
    
    
    
    test = read(readers{1},1);
    test = test(:,:,1);
    s = size(test);
    
    firstFrames = zeros(2,1);
    if mod(s(1),2) == 1
        s(1) = s(1) - 1;
        firstFrames(1) = 2;
    else
        firstFrames(1) = 1;
    end
    
    if mod(s(2),2) == 1
        s(2) = s(2) - 1;
        firstFrames(2) = 2;
    else
        firstFrames(2) = 1;
    end
    
    %keyboard
    % added try/catch
    images = cell(N,1);
    fprintf(1,'\t Loading Images\n');

    for i=1:N
        images{i} = read(readers{groups(i,1)},[groups(i,2) groups(i,3)]);
        images{i} = squeeze(images{i}(:,:,1,:));

    end
    
        
    
    
    
    
    fprintf(1,'\t Making Images\n');    
    parfor i=1:runLength
    
        Z = uint8(zeros(s(1)*subY,s(2)*subX)+255);
   
        for j=1:N
        
            yval = floor((j-1)/subX);
            xval = j - yval*subX - 1;
            
            q = images{j}(firstFrames(1):end,firstFrames(2):end,mod(i-1,lengths(j))+1);
            ss = size(q);
            if ss(1) < s(1) || ss(2) < s(2)
                qq = zeros(s);
                qq(1:ss(1),1:ss(2)) = q;
                q = qq;
            end
            Z((1:s(1)) + yval*s(1),(1:s(2)) + xval*s(2)) = q;
                        
        end
        
        for j=(N+1):(subX*subY)
            yval = floor((j-1)/subX);
            xval = j - yval*subX - 1;
            Z((1:s(1)) + yval*s(1),(1:s(2)) + xval*s(2)) = 0;
        end
            
        a = num2str(i);
        if length(a) < numZeros
            a = [repmat('0',1,numZeros - length(a)) a];
        end
        
        if complementImage
            Z = imcomplement(Z);
        end
        
        
        imwrite(Z,[file_path a '.tif'],'tiff');
        
      
    end

    
    clear images
    