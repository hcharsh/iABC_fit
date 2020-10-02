function [PostDist] = posterior(sort_PRM, lb, ub, n_bins)
    thres = size(sort_PRM, 1); % number of best parameter combinations to be selected
    num_prm = size(sort_PRM, 2); % number of free parameters
    
    % output: binned frequency function of the parameter values selected
    PostDist = zeros(n_bins, num_prm);
    
    for ind_iter = 1:thres
        clear ind_interest;
        for ind_prm = 1:num_prm
            clear abs_rlz rlz;
            
            % current parameter value
            abs_rlz = sort_PRM(ind_iter, ind_prm);
            
            % abs_rlz normalized according to lb->0 and ub->1
            rlz = (abs_rlz - lb(ind_prm, 1))/(ub(ind_prm, 1) - lb(ind_prm, 1));
            
            % bin index for rlz obtained
            bin_rlz = floor(n_bins*rlz) + 1;
            
            % frequency function updated
            PostDist(bin_rlz, ind_prm) = PostDist(bin_rlz, ind_prm) + (1/thres);
        end
    end
end

