% opt_t2c_sim_sd2s
%
% this script shows an example of running many simulations with the
% opt_t2c_sim function to assess the effect of systematically changing a 
% parameter controlling metacognitive sensitivity on the simulation results. 
% in this case, we use the included function type2noisySignalLoss as the
% type 2 model and systematically vary the parameter sd2, which controls
% the standard deviation of type 2 noise. 
%
% this scripts mirrors the analaysis of the type 2 noise model conducted in 
% the main manuscript, as displayed in Figure 9.
%
% this script can be readily adapted to run different kinds of simulations. 
% to change details of the simulation, make appropriate edits in the
% section titled "simulation setup". note that by default, settings.ntrials 
% is set to a low value of 1e5 to allow for fast but noisy simulations. to
% increase simulation precision, it may be necessary to substantially
% increase settings.ntrials.
%
% the code contains an option for saving simulated data and/or plots to disk,
% since high-fidelity simulations can have long run times.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

clear

addpath(genpath('..'));

tic

%% simulation setup
% prepare inputs to opt_t2c_sim and set options for saving simulation results

%%% type %%%
% determine simulation type. 
% valid strings are 't2acc', 't2reward', 'calibration', and 'hr2-far2'
type = 'calibration';


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

% t2model.sd2 and t2model.k are defined below, within the loop that iterates over the sd2s variable
% here we define the array of sd2 values over which to iterate
sd2s = 0:.2:2;                % sd2 is std dev of type 2 noise, must be >= 0.
ks   = zeros(size(sd2s));     % k is signal loss, must be in [0,1]. 0 --> no signal loss, 1 --> complete signal loss.

% this t2titlestr is not for use with the plot made by opt_t2c_sim, but
% rather is for use with the plots created at the bottom of this script
t2titlestr = 'using type 2 noise model, varying \sigma_2';


%%% settings %%%
% define simulation settings and options
settings.ntrials       = 1e5;
settings.compute_metad = 0;

% settings.makePlot controls single-simulation plots as output by
% opt_t2c_sim, which we're turning off here or else we'll get length(sd2s)
% single-simulation plots! instead, plots that aggregate over all
% simulations are made later on in this script.
settings.makePlot      = 0;


%%% save options %%%
savedir   = 'results/sd2/';
save_data = 0; % save simulation results to a .mat file in savedir
save_plot = 0; % save plots to .fig and .png files in savedir, and automatically close the Matlab plots


%% run simulation

for i = 1:length(sd2s)

    %%% finish defining t2model for the current simulation
    t2model.sd2 = sd2s(i);
    t2model.k   = ks(i);
    
    % run the simulation
    [opt_sim, opt_calc, M] = opt_t2c_sim(type, param, t2model, settings);    

    % save simulation results in a format convenient for subsequent analyses
    c2_rS1_sim(i)      = opt_sim.c2_rS1;
    c2_rS2_sim(i)      = opt_sim.c2_rS2;
    t2perf_rS1_sim(i)  = opt_sim.t2perf_rS1;
    t2perf_rS2_sim(i)  = opt_sim.t2perf_rS2;
    
    c2_rS1_calc(i)     = opt_calc.c2_rS1;
    c2_rS2_calc(i)     = opt_calc.c2_rS2;
    t2perf_rS1_calc(i) = opt_calc.t2perf_rS1;
    t2perf_rS2_calc(i) = opt_calc.t2perf_rS2;

    md(i)     = M.md;
    md_rS1(i) = M.md_rS1;
    md_rS2(i) = M.md_rS2;
    M_rS1(i)  = M.M_rS1;
    M_rS2(i)  = M.M_rS2;
end


%% save data

if save_data
    
    mkdir(savedir);

    % save overall simulation results
    filename = ['sd2_' type '_d=' num2str(param.d) '_c=' num2str(param.c) '_pS2=' num2str(param.pS2)];
    
    switch type
        case 't2reward'
            Q2 = (param.R.CR2 - param.R.FA2) / (param.R.hit2 - param.R.miss2);
            filename = [filename '_Q2=' num2str(Q2)];
        case 'calibration'
            OT = param.pcorr_thresh / (1 - param.pcorr_thresh);            
            filename = [filename '_OT=' num2str(OT)];
    end
    
    filename = [filename '_nt=1e' num2str(log10(settings.ntrials))];
    save([savedir filename '.mat']);
    
    % save analysis regarding the relationship b/t type 2 sd and M_ratio
    % (which doesn't depend on optimization context)
    if settings.compute_metad
        filename_M = ['sd2_Mratio_d=' num2str(param.d) '_c=' num2str(param.c) '_pS2=' num2str(param.pS2) '_nt=1e' num2str(log10(settings.ntrials))];
        save([savedir filename_M '.mat'], 'param', 'settings', 'md', 'md_rS1', 'md_rS2', 'M_rS1', 'M_rS2');
    end    
end


%% plots

%% plot optimal type 2 criteria as a function of sd2
% plots optimal type 2 criteria as a function of type 2 noise sd2 for each
% response type. the optimal type 2 criterion when ignoring type 2 noise
% and using classical SDT corresponds to the point on the plot where sd2 = 0

yminmax = [min([c2_rS1_sim, c2_rS2_sim]), max([c2_rS1_sim, c2_rS2_sim])];

figure;

% plot "S1" response results
subplot(1,2,1); hold on;
plot(sd2s, c2_rS1_sim, 'o-'); 
ylabel('optimal c_2 for "S1" responses')
xlabel('\sigma_2')
ylim(yminmax)

% plot "S2" response results
subplot(1,2,2); hold on;
plot(sd2s, c2_rS2_sim, 'o-'); 
ylabel('optimal c_2 for "S2" responses')
xlabel('\sigma_2')
ylim(yminmax)

% title
switch type
    case 't2reward'
        Q_2 = (param.R.CR2 - param.R.FA2) / (param.R.hit2 - param.R.miss2);
        titlestr_extra = [', Q_2 = ' num2str(Q_2)];
    case 'calibration'
        OT = param.pcorr_thresh / (1 - param.pcorr_thresh);            
        titlestr_extra = [', O_T = ' num2str(OT)];
    otherwise
        titlestr_extra = [];
end  


titlestr = {['simulation results for ' type ' (n trials = ' num2str(settings.ntrials) ')'], ...
            ['with d'' = ' num2str(param.d) ', c_1 = ' num2str(param.c) ', p(S2) = ' num2str(param.pS2), titlestr_extra], ...
             t2titlestr};
            
try
    sgtitle(titlestr);
catch
    title(titlestr);
end

% save the plot if requested
if save_plot
    saveas(gcf, [savedir filename '_c2.fig'], 'fig');
    saveas(gcf, [savedir filename '_c2.png'], 'png');
    close(gcf);
end


%% plot optimal type 2 performance as a function of sd2
% plots type 2 performance as a function of type 2 noise sd2 for each
% response type, both when taking into account the type 2 model and when 
% ignoring it. when taking the type 2 model into account, optimal type 2
% performance corresponds to the type 2 performance achieved when using the
% simulated optimal type 2 criterion. when ignoring the type 2 model,
% optimal type 2 performance corresponds to the type 2 performance achieved
% when using the type 2 criterion as theoretically derived from classical
% SDT.

yminmax = [min([t2perf_rS1_calc, t2perf_rS2_calc]), max([t2perf_rS1_sim, t2perf_rS2_sim])];

figure; 

% plot "S1" response results
subplot(1,2,1); hold on;
plot(sd2s, t2perf_rS1_sim, 'o-');
plot(sd2s, t2perf_rS1_calc, 'o-');

switch type
    case 't2acc'
        ylabel('optimal type 2 accuracy for "S1" responses')
    case 't2reward'
        ylabel('optimal type 2 reward for "S1" responses')
    case 'calibration'
        ylabel('estimated p(correct) at optimal c_2 for "S1" responses')
    case 'hr2-far2'
        ylabel('optimal HR_2 - FAR_2 for "S1" responses')
end

xlabel('\sigma_2')
legend('accounting for type 2 model (sim)', 'ignoring type 2 model (calc)')
ylim(yminmax)

% plot "S2" response results
subplot(1,2,2); hold on;
plot(sd2s, t2perf_rS2_sim, 'o-'); 
plot(sd2s, t2perf_rS2_calc, 'o-');

switch type
    case 't2acc'
        ylabel('optimal type 2 accuracy for "S2" responses')
    case 't2reward'
        ylabel('optimal type 2 reward for "S2" responses')
    case 'calibration'
        ylabel('estimated p(correct) at optimal c_2 for "S2" responses')
    case 'hr2-far2'
        ylabel('optimal HR_2 - FAR_2 for "S2" responses')
end

xlabel('\sigma_2')
ylim(yminmax)

% title
switch type
    case 't2reward'
        Q_2 = (param.R.CR2 - param.R.FA2) / (param.R.hit2 - param.R.miss2);
        titlestr_extra = [', Q_2 = ' num2str(Q_2)];
    case 'calibration'
        OT = param.pcorr_thresh / (1 - param.pcorr_thresh);            
        titlestr_extra = [', O_T = ' num2str(OT)];
    otherwise
        titlestr_extra = [];
end  

titlestr = {['simulation results for ' type ' (n trials = ' num2str(settings.ntrials) ')'], ...
            ['with d'' = ' num2str(param.d) ', c_1 = ' num2str(param.c) ', p(S2) = ' num2str(param.pS2), titlestr_extra], ...
             t2titlestr};
         
try
    sgtitle(titlestr);
catch
    title(titlestr);
end

% save the plot if requested
if save_plot
    saveas(gcf, [savedir filename '_t2perf.fig'], 'fig');
    saveas(gcf, [savedir filename '_t2perf.png'], 'png');
    close(gcf);
end


%% M_ratio plots

%% plot optimal type 2 criteria as a function of M_ratio
% plots optimal type 2 criteria as a function of M_ratio for each response 
% type. the optimal type 2 criterion when ignoring type 2 noise and using 
% classical SDT corresponds to the point on the plot where M_ratio = 1.
% 
% M_ratio is controlled by sd2 such that as sd2 increases, M_ratio
% decreases.

if settings.compute_metad

    yminmax = [min([c2_rS1_sim, c2_rS2_sim]), max([c2_rS1_sim, c2_rS2_sim])];
    
    figure;
    
    % plot "S1" response results
    subplot(1,2,1); hold on;
    plot(M_rS1, c2_rS1_sim, 'o-');
    ylabel('optimal c_2 for "S1" responses')
    xlabel('M_{ratio} as controlled by \sigma_2')
    ylim(yminmax)
    
    % plot "S2" response results
    subplot(1,2,2); hold on;
    plot(M_rS1, c2_rS2_sim, 'o-');
    ylabel('optimal c_2 for "S2" responses')
    xlabel('M_{ratio} as controlled by \sigma_2')
    ylim(yminmax)
    
    
    % title
    switch type
        case 't2reward'
            Q_2 = (param.R.CR2 - param.R.FA2) / (param.R.hit2 - param.R.miss2);
            titlestr_extra = [', Q_2 = ' num2str(Q_2)];
        case 'calibration'
            OT = param.pcorr_thresh / (1 - param.pcorr_thresh);            
            titlestr_extra = [', O_T = ' num2str(OT)];
        otherwise
            titlestr_extra = [];
    end  
    
    titlestr = {['simulation results for ' type ' (n trials = ' num2str(settings.ntrials) ')'], ...
                ['with d'' = ' num2str(param.d) ', c_1 = ' num2str(param.c) ', p(S2) = ' num2str(param.pS2), titlestr_extra], ...
                 t2titlestr};
                
    try
        sgtitle(titlestr);
    catch
        title(titlestr);
    end
    
    % save the plot if requested
    if save_plot
        saveas(gcf, [savedir filename '_c2_M.fig'], 'fig');
        saveas(gcf, [savedir filename '_c2_M.png'], 'png');
        close(gcf);
    end

end


%% plot optimal type 2 performance as a function of M_ratio
% plots type 2 performance as a function of M_ratio for each response type, 
% both when taking into account the type 2 model and when ignoring it. when 
% taking the type 2 model into account, optimal type 2 performance corresponds 
% to the type 2 performance achieved when using the simulated optimal type 2 
% criterion. when ignoring the type 2 model, optimal type 2 performance 
% corresponds to the type 2 performance achieved when using the type 2 
% criterion as theoretically derived from classical SDT.
%
% M_ratio is controlled by sd2 such that as sd2 increases, M_ratio
% decreases.

if settings.compute_metad

    yminmax = [min([t2perf_rS1_calc, t2perf_rS2_calc]), max([t2perf_rS1_sim, t2perf_rS2_sim])];
    
    figure;
    
    % plot "S1" response results
    subplot(1,2,1); hold on;
    plot(M_rS1, t2perf_rS1_sim, 'o-');
    plot(M_rS1, t2perf_rS1_calc, 'o-');

    switch type
        case 't2acc'
            ylabel('optimal type 2 accuracy for "S1" responses')
        case 't2reward'
            ylabel('optimal type 2 reward for "S1" responses')
        case 'calibration'
            ylabel('estimated p(correct) at optimal c_2 for "S1" responses')
        case 'hr2-far2'
            ylabel('optimal HR_2 - FAR_2 for "S1" responses')
    end    
    
    xlabel('M_{ratio} as controlled by \sigma_2')
    legend('accounting for type 2 model (sim)', 'ignoring type 2 model (calc)')
    ylim(yminmax)
    
    % plot "S2" response results
    subplot(1,2,2); hold on;
    plot(M_rS2, t2perf_rS2_sim, 'o-');
    plot(M_rS2, t2perf_rS2_calc, 'o-');

    switch type
        case 't2acc'
            ylabel('optimal type 2 accuracy for "S2" responses')
        case 't2reward'
            ylabel('optimal type 2 reward for "S2" responses')
        case 'calibration'
            ylabel('estimated p(correct) at optimal c_2 for "S2" responses')
        case 'hr2-far2'
            ylabel('optimal HR_2 - FAR_2 for "S2" responses')
    end    
    
    xlabel('M_{ratio} as controlled by \sigma_2')
    ylim(yminmax)
    
    % title
    switch type
        case 't2reward'
            Q_2 = (param.R.CR2 - param.R.FA2) / (param.R.hit2 - param.R.miss2);
            titlestr_extra = [', Q_2 = ' num2str(Q_2)];
        case 'calibration'
            OT = param.pcorr_thresh / (1 - param.pcorr_thresh);            
            titlestr_extra = [', O_T = ' num2str(OT)];
        otherwise
            titlestr_extra = [];
    end  
    
    titlestr = {['simulation results for ' type ' (n trials = ' num2str(settings.ntrials) ')'], ...
                ['with d'' = ' num2str(param.d) ', c_1 = ' num2str(param.c) ', p(S2) = ' num2str(param.pS2), titlestr_extra], ...
                 t2titlestr};    
    
    try
        sgtitle(titlestr);
    catch
        title(titlestr);
    end
    
    % save the plot if requested
    if save_plot
        saveas(gcf, [savedir filename '_t2perf_M.fig'], 'fig');
        saveas(gcf, [savedir filename '_t2perf_M.png'], 'png');
        close(gcf);
    end

end


%% plot M_ratio as a function of sd2
% plots M_ratio as a function of sd2 for each response type

if settings.compute_metad
    yminmax = [0, 1.1];

    figure;
    
    % plot "S1" response results
    subplot(1,2,1); hold on;
    plot(sd2s, M_rS1, 'o-'); 
    plot(sd2s, ones(size(sd2s)), 'k--');
    ylabel('M_{ratio} for "S1" responses')
    xlabel('\sigma_2')
    ylim(yminmax)

    % plot "S2" response results
    subplot(1,2,2); hold on;
    plot(sd2s, M_rS2, 'o-'); 
    plot(sd2s, ones(size(sd2s)), 'k--');
    ylabel('M_{ratio} for "S2" responses')
    xlabel('\sigma_2')
    ylim(yminmax)


    % title
    titlestr = {['simulation results (n trials = ' num2str(settings.ntrials) ')'], ...
                ['with d'' = ' num2str(param.d) ', c_1 = ' num2str(param.c) ', p(S2) = ' num2str(param.pS2)], ...
                 t2titlestr};               
             
    try
        sgtitle(titlestr);
    catch
        title(titlestr);
    end
    
    % save the plot if requested
    if save_plot
        saveas(gcf, [savedir filename_M '.fig'], 'fig');
        saveas(gcf, [savedir filename_M '.png'], 'png');
        close(gcf);
    end    
end

toc