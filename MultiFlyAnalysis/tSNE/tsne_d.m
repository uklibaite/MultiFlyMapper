function [ydata,betas,P,costs] = tsne_d(D, no_dims, perplexity, relTol, parameters)
%TSNE_D Performs symmetric t-SNE on the pairwise Euclidean distance matrix D
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego
%
% Modifications by Gordon J. Berman, 2014
% Princeton University


    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
    sigmaTolerance = parameters.sigmaTolerance;
    maxNeighbors = parameters.maxNeighbors;
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
    else
        initial_solution = false;
    end
    
    % Compute joint probabilities
    D = D / max(D(:));    
    
    if ischar(perplexity)
        perplexity = str2double(perplexity);
    end
    
    if ischar(relTol)
        relTol = str2double(relTol);
    end
    
        % normalize distances to avoid machine precision errors
    [P,betas] = d2p_sparse(D .^ 2, perplexity, sigmaTolerance,maxNeighbors);    
        % compute affinities using fixed perplexity
    
    
    % Run t-SNE
    if initial_solution
        [ydata,costs] = tsne_p_sparse(P, ydata, relTol, parameters);
    else
        [ydata,costs] = tsne_p_sparse(P, no_dims, relTol, parameters);
    end
    