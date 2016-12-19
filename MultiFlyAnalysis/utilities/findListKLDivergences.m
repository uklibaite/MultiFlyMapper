function [D,entropies] = findListKLDivergences(data,data2)

    %N = length(data(:,1));
    M = length(data2(:,1));
    logData = log(data);
    
    
    entropies = -sum(data.*logData,2);
    clear logData;  

    logData2 = log(data2);  

    D = - data * logData2';
    
    D = bsxfun(@minus,D,entropies); 
    %D = D - repmat(entropies,1,M);
        
    D = D ./ log(2);