% opt_t2c_sim_run
%
% this script shows an example of running a single simulation with the
% opt_t2c_sim function, using the included function type2noisySignalLoss as 
% the type 2 model. for example scripts running many simulations over many
% settings of type 2 model parameters, see opt_t2c_sim_sd2s and
% opt_t2c_sim_ks.
%
% this script can be readily adapted to run different kinds of simulations.
% to change details of the simulation, make appropriate edits in the
% section titled "simulation setup". note that by default, settings.ntrials 
% is set to a low value of 1e5 to allow for fast but noisy simulations. to
% increase simulation precision, it may be necessary to substantially
% increase settings.ntrials.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters


clear

addpath(genpath('..'));


%% simulation setup
% prepare inputs to opt_t2c_sim

%%% type %%%
% determine type simulation type. 
% valid strings are 't2acc', 't2reward', 'calibration', and 'hr2-far2'
type = 't2reward';

%%% param %%%
% define simulation parameters
param.d   = 2;
param.c   = 0;
param.pS2 = 0.5;

% type 2 reward matrix (for when type is 't2reward'; ignored otherwise)
param.R.hit2  = 1;
param.R.miss2 = 0;
param.R.CR2   = 3;
param.R.FA2   = 0;

% p(correct) threshold (for when type is 'calibration'; ignored otherwise)
param.pcorr_thresh = 0.75;


%%% t2model %%%
% define the type 2 model and parameters
t2model.t2fn = @type2noisySignalLoss; % funtion used to define the type 2 model
t2model.sd2  = 0;                     % sd2 is std dev of type 2 noise, must be >= 0.
t2model.k    = 0;                     % k is signal loss, must be in [0,1]. 0 --> no signal loss, 1 --> complete signal loss.

% define the title string used to describe the type 2 model in the plot made by opt_t2c_sim
t2model.t2titlestr = ['using type 2 noisy signal loss model (\sigma_2 = ' num2str(t2model.sd2), ', k = ' num2str(t2model.k) ')'];


%%% settings %%%
% define simulation settings and options
settings.ntrials       = 1e5;
settings.compute_metad = 0;
settings.makePlot      = 1;


%% run the simulation

[opt_sim, opt_calc, M] = opt_t2c_sim(type, param, t2model, settings);
