function PARAMETER = prm_rlz(currDist, lb, ub, n_rlz)
    [n_bins, n_prm] = size(currDist);
    unifDist = linspace(0, 1, n_bins)';
    
    % uniform random samplig using Latin Hyper-cube algorithm
    X1 = lhsdesign(n_rlz, n_prm,'iterations', 5);
    
    % using X1 to sample according to the distribution of interest: currDist
    for ind_prm = 1:n_prm
        % scaling
        low = lb(ind_prm, 1); high = ub(ind_prm, 1); dif = high - low;
        for ind_rlz = 1:n_rlz
            
            % random sample between 0 and 1.
            rlz_unif = X1(ind_rlz, ind_prm);
            
            bin_dist = 1;
            while 1
                if rlz_unif<currDist(bin_dist, ind_prm)
                    break;
                else
                    bin_dist = bin_dist + 1; % bin in which rlz_unif lies
                end
            end            
%             disp([rlz_unif, bin_dist]);
            
            % [d0, d1] is the {bin_dist}th bin in currDist
            d1 = currDist(bin_dist, ind_prm);
            d0 = currDist(bin_dist - 1, ind_prm);
            
            % [x0, x1] is the {bin_dist}th bin in uniform distribution
            x1 = unifDist(bin_dist, 1);
            x0 = unifDist(bin_dist - 1, 1);
            
            % % unif -> dist: converting "rlz_unif" to appropriate value 
            % according to distribution
            rlz_dist = x0 + (x1 - x0)*(rlz_unif - d0)/(d1 - d0);
            val = low + rlz_dist*dif;
            
            PARAMETER(ind_rlz, ind_prm) = val; % scaled parameter combination samples
            clear d0 d1 x0 x1 rlz_unif bin_dist rlz_dist val;
        end
        clear low high dif;
    end
end