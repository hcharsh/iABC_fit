## Cellular life cycle of positive sense RNA virus
MATLAB codes to fit the viral life cycle dynamics using iABC algorithm.

### Supplementary to following article:
"Life cycle process dependencies of positive-sense RNA viruses suggest
strategies for inhibiting productive cellular infection" - Chhajer et al.

Current location:
https://www.biorxiv.org/content/10.1101/2020.09.19.304576v2

Cite the article upon use of the codes.

### Context-independent functions

iABC_virus: iABC algorithm to iteratively improve parameter value
distribution estimate.

posterior: converts the parameter occurence into a distribution function

prm_rlz: generating random samples according to any given distribution.
Sampling based on Latin Hyper cube algorithm.

estimate_cmf: summarizes the estimated distributions in terms of
quantiles.

plot_PracId: characterizing practical identifiability (PI) using pairwise
correlation between the corresponding parameter values estimated in the
final iteration of iABC.

### Context-dependent functions

#### Model specific

indexing: Global information file: parameters, variables and virus index
and related information being repeatedly used in various functions

model_tau: the life cycle model

initial guess: input the lower and upper bounds for parameter values

#### Data specific

Error_HCV: simulates the model and evaluates the error (chi_sq) between
prediction and experimental data. Modify this to change the data file
Error_JEV, Error_PV files are available on request

calc_error: used to pass current parameter combinations to Error_HCV

plot_dyn_HCV: simulates the model using the parameter combinations 
selected; plots the dynamics and compares them with data.

### File to run - users need to specify iABC options

run_this: executes the iABC_virus by passing iABC algorithm options
