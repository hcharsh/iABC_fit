%% variable indexing
cRNAind     = 1;
NSPind      = 2;
SPind       = 3;
RCind       = 4;
rcRNAind    = 5;
tot_Vind    = 6;
tot_ind     = 6;

%% parameter indexing
kt_ind      = 1;
krc_ind     = 2;
tau_ind     = 3;

kr_ind      = tau_ind + 1;
kexport_ind = tau_ind + 2;
rcsat_ind   = tau_ind + 3;

ka_ind      = tau_ind + 4;

dgrd_r_ind  = ka_ind + 1;
dgrd_p_ind  = ka_ind + 2;
dgrd_v_ind  = ka_ind + 3;

%% parameter labels
% prm_name - symbol; img_name - name suitable for naming images

prm_name(1, kt_ind) = {'k_t'};
img_name(1, kt_ind) = {'kT'};

prm_name(1, krc_ind) = {'k_{c}'};
img_name(1, krc_ind) = {'kRC'};

prm_name(1, tau_ind) = {'\tau_F'};
img_name(1, tau_ind) = {'tau'};

prm_name(1, kr_ind) = {'k_r'};
img_name(1, kr_ind) = {'kr'};

prm_name(1, kexport_ind) = {'k_{e}'};
img_name(1, kexport_ind) = {'kexport'};

prm_name(1, rcsat_ind) = {'N_{C}'};
img_name(1, rcsat_ind) = {'RCsat'};

prm_name(1, ka_ind) = {'k_a'};
img_name(1, ka_ind) = {'kA'};

%% virus names
virus_name(1, 1) = {'HCV'};
virus_name(1, 2) = {'JEV'};
virus_name(1, 3) = {'PV'};

%% color scheme for variables
col_psr = [1 0 0]; % (+)vRNA in cell
col_tsr = [0.75 0 1]; % total vRNA in cell
col_msr = [0 0 0.75]; % (-)vRNA in cell
col_rdrp_n  = [0 1 0.2]; % normalized RdRp. activity
col_iv  = [0 0 0]; % virus
col_sp = 0.3*[3 2 0]; % str. protein

