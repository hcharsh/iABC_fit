function estimate_cmf(virus_ind)
    
    % virus_ind indicates which result file is to analyzed. The file info
    % is present in file "indexing"; 
    indexing;
    
    % file to be used for analysis
    load(['session_main_', virus_name{1, virus_ind}, '.mat'], 'Pdf_all', 'lb', 'ub');
    % Pdf_all: normalized output of the iABC estimation
    % lb and ub (vectors): the lower and upper bounding the values of each
    % parameter
    
    % file in which the distribution summary is stored
    session_cmf_name = ['cmf_estimate_', virus_name{1, virus_ind}];

    %% distribution info
    % n_bins: number is bins in which the distribution is discretized
    % n_prm: number of free parameters that were estimated
    % n_nest_1: number of iterations in iABC + 1
    [n_bin, n_prm, n_nest_1] = size(Pdf_all);
    LB_info = lb(:, 1); % copying the info
    UB_info = ub(:, 1); % copying the info
    PD_last = Pdf_all(:, :, n_nest_1); % the output of the final iteration

    %% working
    CMF = zeros(n_bin, n_prm);
    CMF(1, :) = PD_last(1, :);
    for bin_iter = 2:n_bin
        CMF(bin_iter, :) = PD_last(bin_iter, :) + CMF(bin_iter - 1, :);
    end
    
    DELTA_PRM = UB_info - LB_info; % for scaling the values
    
    %% find 10%, 25%, 50% (median), 75% and 90% quantiles of the distribution
    f_arr = [0.1; 0.25; 0.5; 0.75; 0.9];
    
    for ind_prm = 1:n_prm
        for ind = 1:size(f_arr, 1)
            f_val = f_arr(ind, 1);
            bin_iter = 1;
            while bin_iter < 1+n_bin
                if CMF(bin_iter, ind_prm) > f_val
                    break;
                end
                bin_iter = bin_iter + 1;
            end
            if bin_iter == 1
                prev_val = 0;
            else
                prev_val = CMF(bin_iter - 1, ind_prm);
            end
            slop = CMF(bin_iter, ind_prm) - prev_val;
            bin_num = (f_val - prev_val)/slop + bin_iter - 1;
            est = LB_info(ind_prm, 1) + (bin_num)*DELTA_PRM(ind_prm, 1)/(n_bin);
            Virus_est(ind_prm, ind) = est;
        end
    end
    save([session_cmf_name, '.mat'], 'Virus_est');
%     close all;
end
