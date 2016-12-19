function [D,entropies] = findKLDivergences(data)

    N = length(data(:,1));
    logData = log(data);
    logData(isinf(logData) | isnan(logData)) = 0;
    
    entropies = -sum(data.*logData,2);
    
    D = - data * logData';
    D = bsxfun(@minus,D,entropies);
    %for i=1:N
    %    D(i,:) = D(i,:) - entropies(i);
    %end
    
    D = D ./ log(2);
    D(1:(N+1):end) = 0;