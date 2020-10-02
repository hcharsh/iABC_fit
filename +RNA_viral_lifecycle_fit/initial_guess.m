function [lb, ub] = initial_guess()

    indexing; % info about parameter, variable and virus indexing
    % lb: vector - lower bound of initial guess for each free parameter
    % ub: vector - upper bound of initial guess for each free parameter
    
    % in this case we have 7 free parameters, and we fit them in the log
    % scale (base 10) in order to explore a large dynamic range
    
    lb(kt_ind, 1) = 0;
    ub(kt_ind, 1) = 3;

    lb(krc_ind, 1) = -5;
    ub(krc_ind, 1) = -1;

    lb(tau_ind, 1) = 0;
    ub(tau_ind, 1) = 2;

    lb(kr_ind, 1) = -1;
    ub(kr_ind, 1) = 3;

    lb(kexport_ind, 1) = -2;
    ub(kexport_ind, 1) = 2;

    lb(rcsat_ind, 1) = 1;
    ub(rcsat_ind, 1) = 4;

    lb(ka_ind, 1) = -10;
    ub(ka_ind, 1) = 0;

end
