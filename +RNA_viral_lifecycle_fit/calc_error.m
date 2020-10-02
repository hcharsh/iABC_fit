function err_vec = calc_error(parameter_matrix, virus_ind)

    % % function inputs
    % parameter_matrix: matrix of sampled parameter combinations
    % virus_ind identifies: (1) HCV; (2) JEV; (3) PV.
    % % function output
    % err_vec: stores error calculated for each sample point
    % % extra function used here:
    % Error_HCV (user defined) - calculates the deviation between model
    % prediction and HCV dynamics
    % Error_JEV, Error_PV files are available on request
    
    % n_rlz: number of parameter combination samples to be analyzed
    n_rlz = size(parameter_matrix, 1);
    
    err_vec = zeros(n_rlz, 1);
    
    for iter = 1:n_rlz
        clear prm_arr prm WSSE;
        
        prm_arr = parameter_matrix(iter, :);
        
        prm = prm_arr';
        
        if mod(iter-1, n_rlz/4)==0
            disp([num2str(100*(iter-1)/n_rlz), '% of current iteration completed']);
            % display completion status
        end
        
        % evaluating error between prediction and data for each parameter
        % combination for (1) HCV, (2) JEV or (3)PV.
        if virus_ind == 1
            WSSE = Error_HCV(prm);
        elseif virus_ind == 2
            WSSE = Error_JEV(prm);
        elseif virus_ind == 3
            WSSE = Error_PV(prm);
        end
        
        err_vec(iter, :) = WSSE;
    end    
end