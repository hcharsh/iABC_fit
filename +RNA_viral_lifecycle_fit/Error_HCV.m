function [WSSE] = Error_HCV(parameters)
% % calculates the deviation between model prediction and HCV dynamics

%% data: taken from
% Aunins et al, J. Virol (2018)

    totRNA_obs = [... time log10(vRNA)
        9	0.935
        12	1.065
        15	1.65
        18	1.88
        21	2.405
        24	3.255
        30	3.555
        36	3.67
        ];

    protein_obs = [... time log10(NSP) log10(SP)
        3	1.41	3.05
        6	2.01	2.78
        9	2.52	2.51
        12	2.62	2.56
        15	2.39	2.22
        18	2.48	2.96
        24	3.04	3.29
        30	3.67	4.16
        36	3.86	4.43
        ];
    
    virus_obs = [... time mean sd   (in V_T levels in linear scale)
        3	3.72	0.76
        6	2.725	0.745
        9	2.54	0.21
        12	2.265	0.135
        15	1.79	0.77
        18	1.5	0.33
        21	0.865	0.335
        24	1.025	0.555
        30	2.12	0.78
        36	6.91	1.07
        42	12.525	6.495
        48	22.615	5.485
        ];

%% options
    tspan = [0, 48]; % simulation time span
    global nSP; nSP = 180; % number of str. protein per virion for HCV

%% fixed parameters, initial conditions and parameter augentation

    % fixed parameters: degradation rates
    deg_rna = 0.25; % of cyt. RNA
    deg_protein = 0.11; % of SP, NSP in cyt.
    deg_virus = 0.063; % of V_T
    
    % initial conditions based on "Aunins et al, J. Virol (2018)"
    MOI_virus = 4.6; % for V_T
    MOI_rna = 3; % for RNA in cyt.
    
    % parameters augmentation
    prm = [parameters; log10(deg_rna); log10(deg_protein); log10(deg_virus)];
    
%% data: time, mean of observation and std. deviation
    time_totRNA = totRNA_obs(:, 1);
    totRNA_obs_log_mean = totRNA_obs(:, 2);
    totRNA_obs_log_sd = log10(1.25); % assuing 25% relative error for RNA level
    
    time_inf_virus = virus_obs(:, 1);
    inf_virus_obs_linear_mean = virus_obs(:, 2);
    inf_virus_obs_linear_sd = virus_obs(:, 3); % assuing reported error
    
    time_SP = protein_obs(:, 1);
    SP_obs_log_mean = protein_obs(:, 3);
    SP_obs_log_sd = log10(1.25); % assuing 25% relative error for str. protein level
        
%% simulation

    indexing;    
    % initial conditions
    y0 = zeros(tot_ind, 1);
    y0(cRNAind) = MOI_rna;
    y0(SPind) = MOI_rna*180;
    y0(tot_Vind) = MOI_virus;
    
    [t_iter,y_iter] = ode23s(@(t_iter,y_iter) model_tau(prm, t_iter, y_iter), tspan, y0);
    
    pred_totRNA = y_iter(:, cRNAind) + y_iter(:, rcRNAind) + 2*y_iter(:, RCind); % both + and - ssRNA
    p_totRNA = interp1(t_iter, pred_totRNA, time_totRNA);
    % comparison of pred. and data in the log scale
    sum_diff_totRNA_log = 0;
    for ind_t = 1:size(time_totRNA, 1)
        diff_totRNA_log = (log10(p_totRNA(ind_t, 1)) - totRNA_obs_log_mean(ind_t, 1));
        sum_diff_totRNA_log = sum_diff_totRNA_log + ((diff_totRNA_log/totRNA_obs_log_sd)^2);
    end


    pred_SP = y_iter(:, SPind);
    p_SP = interp1(t_iter, pred_SP, time_SP);
    % comparison of pred. and data in the log scale
    sum_diff_SP_log = 0;
    for ind_t = 1:size(time_SP, 1)
        diff_SP_log = (log10(p_SP(ind_t, 1)) - SP_obs_log_mean(ind_t, 1));
        sum_diff_SP_log = sum_diff_SP_log + ((diff_SP_log/SP_obs_log_sd)^2);
    end

    pred_V = y_iter(:, tot_Vind);
    p_V = interp1(t_iter, pred_V, time_inf_virus);
    % comparison of pred. and data in the linear scale
    sum_diff_inf_V_linear = 0;
    for ind_t = 1:size(time_inf_virus, 1)
        wt_diff_inf_V_linear = ( p_V(ind_t, 1) - inf_virus_obs_linear_mean(ind_t, 1) )/inf_virus_obs_linear_sd(ind_t, 1);
        sum_diff_inf_V_linear = sum_diff_inf_V_linear + (wt_diff_inf_V_linear^2);
    end
    
    % total error between predicted and observed dynamics
    WSSE = sum_diff_totRNA_log + sum_diff_SP_log + sum_diff_inf_V_linear;
    
end