function [ydata,costs] = tsne_p_sparse(P, no_dims, relTol, parameters)
%TSNE_P Performs symmetric t-SNE on transition matrix P
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
    
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    n = size(P, 1);
    
    momentum = parameters.momentum;
    final_momentum = parameters.final_momentum;
    mom_switch_iter = parameters.mom_switch_iter;
    stop_lying_iter = parameters.stop_lying_iter;
    max_iter = parameters.max_iter;
    epsilon = parameters.epsilon;
    lie_multiplier = parameters.lie_multiplier;
    min_gain = parameters.min_gain;
    intial_variation =relTol; %%%%!@$$
    parameters.initial_variation = relTol;
    %readout = parameters.readout;
    readout=1;

    
    old_cost = 1e10;
    costs = zeros(max_iter,1);
    
    % Make sure P-vals are set properly
    P(1:(n + 1):end) = 0;                      % set diagonal to zero
    P = 0.5 * (P + P');                        % symmetrize P-values
    idx = find(P > 0);
    
    P(idx) = P(idx) ./ sum(P(idx));            % make sure P-values sum to one
    const = sum(P(idx) .* log2(P(idx)));        % constant in KL divergence
    
    if ~initial_solution
        P = P * lie_multiplier;                % lie about the P-vals to find better local minimum
        lying_stopped = false;
    else
        lying_stopped = true;
    end
    
    % Initialize the solution
    if ~initial_solution
        ydata = intial_variation * randn(n, no_dims);
    end
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));

    % Run the iterations
    for iter=1:max_iter
        tic
               
        Q = 1 ./ (1 + findAllDistances(ydata).^2);
        Q(1:n+1:end) = 0;
        Z = sum(Q(:));
        Q = Q./Z;
        
        % Compute the gradients (faster implementation)
        L = Z * (P - Q) .* Q;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
        
        % Update the solution (note that the y_grads are actually -y_grads)
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...    
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));

        cost = const - sum(P(idx) .* log2(Q(idx)));
        costs(iter) = cost;
        diffVal = (old_cost - cost) / old_cost;
        old_cost = cost;
        
        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
            lying_stopped = true;
        end
        
        if iter == stop_lying_iter && ~initial_solution
            P = P ./ lie_multiplier;     
        end
        
        % Print out progress
        if ~rem(iter, readout)
            
            disp(['Iteration ' num2str(iter) ': error is ' ...
                num2str(cost) ,', change is ' num2str(diffVal)]);
        end
        

        if abs(diffVal) < relTol && lying_stopped && iter > 10
            break;
        end
        
    end
    
    
    costs = costs(1:iter);
    
    