function folders = findAllFoldersInFolders(folder,fileType,frontConstraint)


MAXIMAGES = 10000000;

if nargin==1
    fileType = '.tif';
end

if nargin < 3 || isempty(frontConstraint) == 1
    frontConstraint = '';
end


[status,files] = unix(['ls ' folder]);

    

if status == 0
    folders = regexp(regexp(files,'\t','split'),'\n','split');
    f = folders{1}';
    for i=2:length(folders)
        f = [f; folders{i}'];
    end
    folders = sort(f);
    while isempty(folders{1}) == 1
        folders = folders(2:end);
    end
    
else
    folders = {};
end