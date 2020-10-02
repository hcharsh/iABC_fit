function plot_dyn_HCV()
warning off;
virus_ind = 1;
indexing;

%% figure formatting
ms = 15; lw_ind = 0.2; cs = 10; lw_mean = 2; ms_ind = 2; lw_err = 1.5;
frame_x0 = 0.1; frame_y0 = 0.1; frame_x1 = 0.4; frame_y1 = 0.72;
fs = 24;

%% details and options

img_file_name = [virus_name{1, virus_ind}, '_sim'];
% image file name to save the model fit

load(['session_main_', virus_name{1, virus_ind}, '.mat'], 'PRM_nest');
% iABC file used for analysis

% model simulation time span
tspan = [0, 48];

% x-axis limits
xs = tspan(1, 1);
xe = tspan(1, 2);
% y-axis limits
ys = -2; ye = 6;

% dividing the time axis into finer intervals
n_tim = 21;
t_fin = linspace(xs, xe, n_tim);

% number of str. protein per one virus particle
global nSP; nSP = 180;

% % % % data taken from Aunins et al, J. Virol (2018)

% total vRNA per cell
totRNA_obs = [...
    9	0.935
    12	1.065
    15	1.65
    18	1.88
    21	2.405
    24	3.255
    30	3.555
    36	3.67
    ];

% NSP And SP (third col) per cell
protein_obs = [... time NSP SP
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

% virus produced per cell
virus_obs = [... time mean sd
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


% (see Error_HCV file for the undefined notations)

deg_rna = 0.25;
deg_protein = 0.11;
deg_virus = 0.063;
MOI_virus = 4.6;
MOI_rna = 3;


% copying data in new variables

time_totRNA = totRNA_obs(:, 1);
totRNA_obs_log_mean = totRNA_obs(:, 2);
totRNA_obs_log_sd = log10(1.25)*ones(size(time_totRNA));

time_inf_virus = virus_obs(:, 1);
inf_virus_obs_linear_mean = virus_obs(:, 2);
inf_virus_obs_linear_sd = virus_obs(:, 3);
% err_virus_neg, err_virus_pos converts linear std. error into
% approprite log scale
err_virus_neg = log10(inf_virus_obs_linear_mean) - log10(inf_virus_obs_linear_mean - inf_virus_obs_linear_sd);
err_virus_pos = log10(inf_virus_obs_linear_mean + inf_virus_obs_linear_sd) - log10(inf_virus_obs_linear_mean);

time_SP = protein_obs(:, 1);
SP_obs_log_mean = protein_obs(:, 3);
SP_obs_log_sd = log10(1.25)*ones(size(time_SP));


% fitting options and results
thres = size(PRM_nest, 1); % number of selected parameter combinations
n_nest_1 = size(PRM_nest, 3); % number of iterations + 1
PRM_last = PRM_nest(:, :, n_nest_1); % selected parameter combinations

% figure initialization
f = figure('Units', 'normalized', 'Position',[frame_x0 frame_y0 frame_x1 frame_y1]); hold on;


% plotting data
errorbar(time_totRNA, totRNA_obs_log_mean, totRNA_obs_log_sd, 'o',...
    'Color', col_tsr, 'MarkerSize', ms, 'Linewidth', lw_err);

errorbar(time_SP, SP_obs_log_mean, SP_obs_log_sd, 'o',...
    'Color', col_sp, 'MarkerSize', ms, 'Linewidth', lw_err);

errorbar(time_inf_virus, log10(inf_virus_obs_linear_mean), err_virus_neg, err_virus_pos, 'o',...
    'Color', col_iv, 'MarkerSize', ms, 'Linewidth', lw_err);   

% initial condition
y0 = zeros(tot_ind, 1);
y0(cRNAind) = MOI_rna;
y0(SPind) = MOI_rna*180;
y0(tot_Vind) = MOI_virus;

%% simulation
% here we simulate model for each selected parameter combination; we plot
% dynamics predicted from each combination alongwith their average.

% initializing arrays to store average dynamics
PTR_mean = zeros(1, n_tim); % total vRNA
PSP_mean = zeros(1, n_tim); % str. protein
PIV_mean = zeros(1, n_tim); % virus produced

for iter = 1:thres
    clear prm t_iter y_iter pred_ssRNA pred_msRNA pred_Inf_virus psr pmr piv;

    parameters = PRM_last(iter, :)'; % current combination
    
    % parameters augmentation
    prm = [parameters; log10(deg_rna); log10(deg_protein); log10(deg_virus)];

    % simulating the dynamics
    [t_iter,y_iter] = ode23s(@(t_iter,y_iter) model_tau(prm, t_iter, y_iter), tspan, y0);
    
    % extracting variables of interest and interpolating the dynamics for
    % ease of averaging and plotting
    pred_totRNA = y_iter(:, cRNAind) + y_iter(:, rcRNAind) + 2*y_iter(:, RCind); % both + and - ssRNA
    ptr = interp1(t_iter, pred_totRNA, t_fin);
    
    pred_SP = y_iter(:, SPind);
    psp = interp1(t_iter, pred_SP, t_fin);
    
    pred_virus = y_iter(:, tot_Vind);
    piv = interp1(t_iter, pred_virus, t_fin);
    
    % plotting prediction for each selected parameter combination
    plot1 = plot(t_fin, log10(ptr), '-',...
        t_fin, log10(psp), '-',...
        t_fin, log10(piv), '-',...
        'LineWidth', lw_ind,...
        'Markersize', ms_ind);
    
    % assigning color based on varibales
    set(plot1(1), 'Color', 0.7+0.3*col_tsr);
    set(plot1(2), 'Color', 0.7+0.3*col_sp);
    set(plot1(3), 'Color', 0.7+0.3*col_iv);
    
    % updating average dynamics
    PTR_mean = PTR_mean + ptr/thres;
    PSP_mean = PSP_mean + psp/thres;
    PIV_mean = PIV_mean + piv/thres;
    
    clear plot1;
end

hold on;
% plotting the average
plot_mean = plot(t_fin, log10(PTR_mean), '-',...
    t_fin, log10(PSP_mean), '-',...
    t_fin, log10(PIV_mean), '-',...
    'LineWidth', lw_mean);
% assigning color based on varibales
set(plot_mean(1), 'Color', col_tsr);
set(plot_mean(2), 'Color', col_sp);
set(plot_mean(3), 'Color', col_iv);

% re-plotting the data, so that they are clear
errorbar(time_totRNA, totRNA_obs_log_mean, totRNA_obs_log_sd, 'o',...
    'Color', col_tsr, 'MarkerSize', ms, 'Linewidth', lw_err);

errorbar(time_SP, SP_obs_log_mean, SP_obs_log_sd, 'o',...
    'Color', col_sp, 'MarkerSize', ms, 'Linewidth', lw_err);

errorbar(time_inf_virus, log10(inf_virus_obs_linear_mean), err_virus_neg, err_virus_pos, 'o',...
    'Color', col_iv, 'MarkerSize', ms, 'Linewidth', lw_err);   

% figure annotation
legend_names = {'Total vRNA', 'P_{S}', 'Viruses'};
xlabel('Time (h)');
ylabel('Molecules per cell');
title('HCV');%, 'Fontsize', fs-8);
ax1 = gca;
ax1.XTick = linspace(xs, xe, 4);
ax1.YTick = [-2:2:6];
ax1.YTickLabel = {'10^{-2}', '1^{}', '10^{2}', '10^{4}', '10^{6}'};
set(ax1, 'Fontsize', fs);
legend(legend_names, 'Fontsize', fs-4, 'Location', 'NorthWest');
legend boxoff;

% fixing range of axes
sp = xe/20;
xlim([xs-sp xe+sp]);
ylim([ys ye]);

% saving the plot
print(f, img_file_name,'-djpeg','-r960') 
end