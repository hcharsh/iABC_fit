function dydt = model_tau(parameters, t, y)

    % 'parameters' contains the parameter values in log (base 10) scale.
    prm = 10.^parameters;
    
    indexing; % info about parameter, variable and virus indexing
    
    global nSP; % number of str. protein per virion for the virus (defined in Error_HCV)
    
    %% assigning values to parameters
    kt = prm(kt_ind,1);
    krc = prm(krc_ind,1);
    
    tau = prm(tau_ind,1);
    
    kr = prm(kr_ind,1);
    k_export = prm(kexport_ind,1);
    rcsat = prm(rcsat_ind,1);
    
    ka = prm(ka_ind,1);
    
    dgrd_r = prm(dgrd_r_ind,1);
    dgrd_p = prm(dgrd_p_ind,1);
    dgrd_v = prm(dgrd_v_ind,1);

    %% rate evaluations according to model
    
    Protein_prdtion = kt*y(cRNAind); % protein synthesis
    
    RNA_RC2cyt_rate = k_export*y(rcRNAind); % +RNA export from CM to cyt.
    
    f_CM   = 1 - exp(-(t/tau)^4); % normalized dynamics of CM formation
    
    RC_form_rate    = krc * y(cRNAind) * y(NSPind) *( (f_CM) - (  y(RCind)/rcsat  ) );
    % rate of formation os compartmentalized replication complexes
    
    V_form_rate     = ka*y(cRNAind)*y(SPind); % virus generation rate
    
    %% rate of evolution of variables according to model
    
    dydt(cRNAind, 1)    = RNA_RC2cyt_rate - V_form_rate - dgrd_r*y(cRNAind) - RC_form_rate;
    
    dydt(NSPind, 1)     = Protein_prdtion - RC_form_rate - dgrd_p*y(NSPind);
    dydt(SPind, 1)      = Protein_prdtion - ka*nSP*y(cRNAind)*y(SPind) - dgrd_p*y(SPind);

    dydt(RCind, 1)      = RC_form_rate;
    
    dydt(rcRNAind, 1) = kr*y(RCind)  - RNA_RC2cyt_rate;

    dydt(tot_Vind, 1) = V_form_rate - dgrd_v*y(tot_Vind);


end