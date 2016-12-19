function idx = makeManyBradyMovies_avi(file_path,movies,groups,subX,subY,firstGroup,complementImage)


    
    if nargin < 4 || isempty(subX)
        subX = 8;
    end
    
    if nargin < 5 || isempty(subY)
        subY = 8;
    end
    
    if nargin < 6 || isempty(firstGroup)
        firstGroup = 1;
    end
    
    if nargin < 7 || isempty(complementImage)
        complementImage = true;
    end
    
    
    [status,~]=unix(['ls ' file_path]);
    if status == 1
        unix(['mkdir ' file_path]);
    end
    
    L = length(groups);
    nDigits = max([1,ceil(log(L+1)./log(10))]);
    idx = cell(L,1);
    
    fprintf(1,'Initializing Readers\n');
    N = length(movies);
    readers = cell(N,1);
    parfor i=1:N
        i
        readers{i} = VideoReader(movies{i});
    end
    
    
    for i=firstGroup:L
        
        fprintf(1,'Making Movie #%2i\n',i);
        
        if ~isempty(groups{i})
            
            q = num2str(i);
            new_folder = [file_path repmat('0',1,nDigits - length(q)) q '/'];
            
            [status,~]=unix(['ls ' new_folder]);
            if status == 1
                unix(['mkdir ' new_folder]);
            else
                unix(['rm ' new_folder '*.tif']);
            end
            
            %a = randperm(Z,numImages);
            
            %idx{i} = a;
            
            idx{i} = makeBradyMovie_avis_readers(new_folder,readers,groups{i},subX,subY,complementImage);
            
        end
        
    end