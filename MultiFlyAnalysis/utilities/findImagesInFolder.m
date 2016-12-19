function images = findImagesInFolder(folder,fileType,frontConstraint)

MAXIMAGES = 2000000;

if nargin == 1
    fileType = '.tif';
end

if nargin < 3 || isempty(frontConstraint) == 1
    frontConstraint = '';
end

%[status,files]=unix(['ls ' folder '*' fileType]);

%['ls ' folder '*' fileType]
%status

if nargin < 3 || isempty(frontConstraint) == 1
    [status,files] = unix(['find ' folder ' -name "*' fileType '"']);
else
    [status,files] = unix(['find ' folder ' -name "' frontConstraint '*' fileType '"']);
end


if status == 0
    
    filesOut = regexp(regexp(files,'\t','split'),'\n','split');
    images = cell(MAXIMAGES,1);
    count = 1;
    %images = filesOut{1}';
    for i=1:length(filesOut)
        images(count:(count+length(filesOut{i})-1)) = filesOut{i}';
        count = count + length(filesOut{i});
    end
    images = images(1:(count-1));
    images = sort(images);
    
    if max(returnCellLengths(images)) > 0
        
        while isempty(images{1}) == 1
            images = images(2:end);
        end
    else
        images = [];
    end
    
else
    images = [];
end