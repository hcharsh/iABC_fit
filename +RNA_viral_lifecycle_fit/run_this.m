%% clearing the workspace
clear all; clc; close all; warning off; format shortG;

%% options for the simulations

% number of iterations the estimation refinement occurs: k (SI SM1)
n_nest = 8;

% number of parameters sampled at every iteration: V (SI SM1)
n_rlz = 1e4;

% number of parameter combinations with lowest chi_sq values used to
% estimate the parameter value distribution at each iteration: M (SI SM1)
cut_off = max(floor(0.025*n_rlz), 25);

% number of bins in which the parameter value range is divided to define
% the distribution
n_bins = 50;

% strength of uniform noise to be added: to increase the robustness of
% fitting algorithm
alpha = 0.25;

% Global information file: parameters, variables and virus index and label
% related information being repeatedly used in various functions
indexing;

%% fitting and analysis

virus_ind = 1;
% corresponds to HCV infection dynamics in Huh7 cells (data from Aunins et
% al., 2018, J. Virol.

% displaying the current analysis
disp(['Fitting: ', virus_name{1, virus_ind}]);

% fitting engine
iABC_virus(virus_ind, n_rlz, cut_off, n_nest, n_bins, alpha);

% displaying the correlation in parameter estimate
plot_PracId(virus_ind);

% display model fit
plot_dyn_HCV();

% summarizing the distribution using quantiles
estimate_cmf(virus_ind);
