function iABC_virus(virus_ind, n_rlz, cut_off, n_nest, n_bins, alpha)

%% information on model
indexing;

%% information on virus and inital guess of the range of parameters for virus
[lb, ub] = initial_guess();
n_prm = size(lb, 1); % number of free parameters
for ind_prm = 1:n_prm
    PriorDist(:, ind_prm) = linspace(0, 1, n_bins + 1)';
    % uniform distribtion between 0 and 1.
    % initial guess is a uniform distribution
end
% currDist: stores the current distribution function; the range of values
% here is scaled between 0 and 1.
currDist = PriorDist;

% prm_rlz(currDist, lb, ub, n_rlz): is used to rescale the range from "0 to
% 1" to "lb and ub", and then sample n_rlz values from this re-scaled
% distribution. Sampling done using Latin Hyper cube sampling


for ind_bin = 1:n_bins
    Pdf_all(ind_bin, :, 1) = PriorDist(ind_bin + 1, :) - PriorDist(ind_bin, :);
    % storing the initial guess: uniform probability mass function between 0 and 1.
end

% rv_vec = zeros(n_bins, n_prm);
% for ind_prm = 1:n_prm
%     clear arr;
%     arr = linspace(lb(ind_prm, 1), ub(ind_prm, 1), n_bins)';
%     rv_vec(:, ind_prm) = arr;
% end

%% distribution biased by error/deviation from observations
for nest_ind = 1:n_nest
%     toc; tic; % to show time used during the iterations
    disp(['In progress, iteration (or nest) #', num2str(nest_ind)]); % showing progress

    clear PostDist;
    % PostDist: temporary variable storing the distribution estimated in
    % this iteration
    
    clear PARAMETER err_vec;
    PARAMETER = prm_rlz(currDist, lb, ub, n_rlz);
    % PARAMETER is a "n_prm * n_rlz" matrix, storing n_rlz samplings of
    % parameter combinations. One combination has value for the n_prm free
    % parameters. The sampling is based on CurrDist rescaled according to
    % "lb and ub"
    
    err_vec = calc_error(PARAMETER, virus_ind);
    % err_vec: stores the deviation between prediction (according to
    % PARAMETER) and data. Data is indirectly mediated by virus_ind and
    % "information in calc_error" and virus_ind. 
    
    clear sort_err index sort_PRM;
    
    % sorting the err_vec; arranging the PARAMETER matrix according to the
    % sorted "index" and selecting "cut_off" many parameter combinations
    % based on lowest deviation from data.
    [Sort_err, index] = sort(err_vec);
    sort_P = PARAMETER(index, :);
    sort_PRM = sort_P(1:cut_off, :);
    
    % storing the first and the "cut_off"^{th} lowest fitting error, so as
    % to see convergence
    ERROR_CUTS(nest_ind, :) = [Sort_err(1, 1), Sort_err(cut_off, 1)];
    
    % storing the "cut_off" best parameter combinations at each iteration
    PRM_nest(:, :, nest_ind) = sort_PRM;

    % storing the frequency function based on the "cut_off" best parameter
    % combinations at each iteration. posterior is used estimate
    % distribution from the chosen parameter values
    clear Post_pdf;
    Post_pdf = posterior(sort_PRM, lb, ub, n_bins); % temporary for each iteration
    Pdf_all(:, :, nest_ind + 1) = Post_pdf; % storing Post_pdf at each iteration
    
%     disp([size(Post_pdf), size(PriorDist)]);
    
    %  converting the frequency function into a distribution function
    PostDist = zeros(n_bins + 1, n_prm);
    for ind_bin = 2:(n_bins + 1)
        PostDist(ind_bin, :) = PostDist(ind_bin - 1, :) + Post_pdf(ind_bin - 1, :);
    end
    % temporary save at each iteration
    save(['session_main_', virus_name{1, virus_ind}, '.mat'])

    % introducing uniform noise so as to make the algorithm more robust
    clear alpha_curr;
    alpha_curr = alpha/((nest_ind+1)); % strength of uniform noise to be added
    currDist = (1 - alpha_curr)*PostDist + alpha_curr*PriorDist;
%     disp(alpha_curr);
end

% saving the session
save(['session_main_', virus_name{1, virus_ind}, '.mat'],...
    'n_rlz', 'ERROR_CUTS', 'alpha', 'PRM_nest', 'Pdf_all', 'lb', 'ub');
end
