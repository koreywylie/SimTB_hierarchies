%%%%%%%%%%%%%%%%%%%%%%%%
% Create simulated as hierarchy of networks,
%  using scripts from SimTB (w/ minor modifications)
%
% Requires:
%  get_fractal_SM.m in local directory, including modified fns.:
%    get_fratal_SM()  generates nested hierarchy of spatial maps
%    simtb_makeSM()   generates spatial maps, modified from simtb fn.
%    simtb_createmask() generates inclusive spatial map encompassing all voxels
%
% Example use, creates then runs simulation
% >> sP = simtb_create_sP('experiment_hierarchy_SM');
% >> simtb_main(sP)
%


%% OUTPUT PARAMETERS
%-------------------------------------------------------------------------------
% Directory to save simulation parameters and output
out_path = pwd;
% Prefix for saving output
prefix = 'SimHier-';
% FLAG to write data in NIFTI format rather than matlab
saveNII_FLAG = 1;
% Option to display output and create figures throughout the simulations
verbose_display = 1;
%-------------------------------------------------------------------------------

%% RANDOMIZATION
%-------------------------------------------------------------------------------
seed = round(sum(100*clock))  % randomizes parameter values
% seed = 3571;                    % choose seed for repeatable simulation
simtb_rand_seed(seed);          % set the seed 
%-------------------------------------------------------------------------------

%% SIMULATION DIMENSIONS
%-------------------------------------------------------------------------------
M  = 100;   % number of subjects    
% nC is the number of components defined below, based on number of levels of hierarchy;
nV = 32;  % 32^2 = 1024 total voxels
nT = 300; % number of time points           
TR = 2;   % repetition time 
%-------------------------------------------------------------------------------



%% Hierarchy Parameters
%-------------------------------------------------------------------------------
% Requires modified scripts simtb_create_sP.m & get_fractal_SM.m, to evaluate params. & create spatial maps
HIER_levels = 3;  % number of levels in hierarchy
HIER_degree = 3;  % size of parcellation at each level, defaults to 2 in get_fractal_SM()
HIER_equal_size = true; % if false, partitions drawn from Dirichlet r.v. w/ param a = (HIER_degree,..,HIER_degree)HIER_groupSM_file = fullfile(out_path, [prefix, '_groupSM.mat']);  % path to file w/ shared group spatial map

nC = sum(HIER_degree .^(1:HIER_levels));  % num. of comps. from above params.
%-------------------------------------------------------------------------------




%% SPATIAL SOURCES
%-------------------------------------------------------------------------------
SM_source_ID = 99 * ones(1,nC);  % param. checked but unused, 31:99 are placeholder GM sources defined in modified simtb_SMsource.m
%-------------------------------------------------------------------------------

%% COMPONENT PRESENCE
%-------------------------------------------------------------------------------
% [M x nC] matrix for component presence: 1 if included, 0 otherwise
SM_present = ones(M, nC);
%-------------------------------------------------------------------------------

%% SPATIAL VARIABILITY
%-------------------------------------------------------------------------------           
% Variability related to differences in spatial location and shape.
SM_translate_x = zeros(M,nC);  % no translation in x
SM_translate_y = zeros(M,nC);  % no translation in y
SM_theta       = zeros(M,nC);  % no rotation
%                Note that each 'activation blob' is rotated independently.
SM_spread = ones(M,nC); % Spread < 1 is contraction, spread > 1 is expansion.
%-------------------------------------------------------------------------------


%% TC GENERATION
%-------------------------------------------------------------------------------
% Choose the model for TC generation.  To see defined models:
% >> simtb_countTCmodels

TC_source_type = ones(1,nC);    % convolution with HRF for most components
% % to make statistical moments of data look more like real data
% TC_source_type([comp_CSF1 comp_CSF2]) = 3; % spike model for CSF

TC_source_params = cell(M,nC);  % initialize the cell structure
% Use the same HRF for all subjects and relevant components
P(1) = 6;    % delay of response (relative to onset)
P(2) = 16;   % delay of undershoot (relative to onset)
P(3) = 1;    % dispersion of response
P(4) = 1;    % dispersion of undershoot
P(5) = 6;    % ratio of response to undershoot
P(6) = 0;    % onset (seconds)
P(7) = 32;   % length of kernel (seconds)
[TC_source_params{:}] = deal(P);
%-------------------------------------------------------------------------------

%% UNIQUE EVENTS
%-------------------------------------------------------------------------------
TC_unique_FLAG = 1; % 1 = include unique events
TC_unique_prob = 0.2*ones(1,nC); % [1 x nC] prob of unique event at each TR, w/o event prob. heterogeneity
% 7/28/2022 --kw-- NOTE: TC_unique_prob_ must be [1 x nC], specification below & in parts of manual incorrect
% TC_unique_prob = repmat(unifrnd(0.2, 0.6, 1, nC),M,1); % [M x nC] prob of unique event at each TR, w/ event heterogeneity
TC_unique_amp  = ones(M,nC) + 0.3 * repmat(randn(1,nC),M,1); % [M x nC] matrix of amplitude of unique events, w/ added heterogeneity 
%-------------------------------------------------------------------------------

%% DATASET BASELINE                                 
%-------------------------------------------------------------------------------
% [1 x M] vector of baseline signal intensity for each subject
D_baseline = 800*ones(1,M); % [1 x M] vector of baseline signal intensity
%-------------------------------------------------------------------------------

%% TISSUE TYPES
%-------------------------------------------------------------------------------
% FLAG to include different tissue types (distinct baselines in the data)
D_TT_FLAG = 0;                    % if 0, baseline intensity is constant 
% D_TT_level = [1.15, 0.8, 1, 1.2]; % TT fractional intensities
% To see/modify definitions for tissue profiles:
% >> edit simtb_SMsource.m
%-------------------------------------------------------------------------------

%% PEAK-TO-PEAK PERCENT SIGNAL CHANGE 
%-------------------------------------------------------------------------------
D_pSC = 3 + 0.25*randn(M,nC);   % [M x nC] matrix of percent signal changes 
%-------------------------------------------------------------------------------

%% NOISE
%-------------------------------------------------------------------------------
D_noise_FLAG = 1;               % FLAG to add rician noise to the data
% [1 x M] vector of contrast-to-noise ratio for each subject
% CNR is distributed as uniform between 0.65 and 2.0 across subjects.  
minCNR = 0.65;  maxCNR = 2;
D_CNR = rand(1,M)*(maxCNR-minCNR) + minCNR; 
%-------------------------------------------------------------------------------
