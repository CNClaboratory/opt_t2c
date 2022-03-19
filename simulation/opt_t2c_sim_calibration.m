% opt_t2c_sim_calibration
%
% this script investigates the type 2 criteria that optimize calibration 
% by plotting the effect of d', c, and p(correct)_T on the optimal criteria 
% under classical SDT as well as the type 2 noisy signal loss model 
% (sd2 > 0 and/or k > 0).
%
% results demonstrate that although sd2 affects the optimal type 2
% criteria, the qualitative patterns whereby optimal type 2 criteria become 
% more liberal with increasing d' and more conservative with increasing
% p(correct)_T still hold.
%
% this script complements the analaysis of Figure 6 in the main manuscript
% by adding consideration of how the plots change under type 2 noise and/or 
% signal loss.
%
% the details of this analysis can be readily changed by appropriately 
% editing the "simulation setup" section. note that by default, 
% settings.ntrials is set to a low value of 1e6 to allow for fast but noisy 
% simulations. to increase simulation precision, it may be necessary to 
% substantially increase settings.ntrials.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

clear

addpath(genpath('..'));

%% simulation setup

%%% param %%%
% define simulation parameters
pS2  = 0.5;         % p(S2)
pcts = [.6, .9];    % p(correct) threshold
cs   = [-1, 0, 1];  % type 1 criteria
ds   = .5 : .5 : 3; % d' values


%%% t2model %%%
% define the type 2 model and parameters
t2model.t2fn = @type2noisySignalLoss; % funtion used to define the type 2 model
t2model.sd2  = 1;                     % std dev of type 2 noise, must be >= 0
t2model.k    = 0;                     % signal loss, must be in [0,1]. 0 --> no signal loss, 1 --> complete signal loss.


%%% settings %%%
% define simulation settings and options
settings.ntrials       = 1e6;

% this script is not designed to conduct meta-d' analysis or to work with
% plots made by opt_t2c_sim
settings.compute_metad = 0;
settings.makePlot      = 0;


%% construct the plot

figure;

i_plot = 0;
for i_pct = 1:2
    for i_c = 1:3
        
        c   = cs(i_c);
        pct = pcts(i_pct);
        OT  = pct / (1-pct);
        
        % compute the optimal type 2 criteria
        for i_d = 1:length(ds)
            d = ds(i_d);
            
            param.d   = d;
            param.c   = c;
            param.pS2 = pS2;
            param.pcorr_thresh = pct;
                        
            [opt_sim, opt_calc, M] = opt_t2c_sim('calibration', param, t2model, settings);
            
            c2_rS1_sim(i_d)  = opt_sim.c2_rS1;
            c2_rS2_sim(i_d)  = opt_sim.c2_rS2;

            c2_rS1_calc(i_d) = opt_calc.c2_rS1;
            c2_rS2_calc(i_d) = opt_calc.c2_rS2;
            
        end
        
        
        i_plot = i_plot + 1;
        
        subplot(2, 3, i_plot); hold on;
        plot(ds, c2_rS2_sim, 'b-', 'LineWidth', 1);
        plot(ds, c2_rS2_calc, 'r-', 'LineWidth', 1);
        plot(ds, cs(i_c)*ones(size(ds)), 'k-')
        plot(ds, c2_rS1_calc, 'r--', 'LineWidth', 1);
        plot(ds, c2_rS1_sim, 'b--', 'LineWidth', 1);
        
        plot(ds, zeros(size(ds)), 'k--')
        
        xlabel('d''')
        xlim([min(ds), max(ds)])
        ylim([-3, 3])
        
        switch i_plot
            case 1
                title(['c_1 = ' num2str(c)]);
                ylabel(['p(correct_1)_T = ' num2str(pct)]);
            case 2
                title(['c_1 = ' num2str(c)]);
                legend(['c^{C*}_{2,"S2"} (\sigma_2=' num2str(t2model.sd2) ', k=' num2str(t2model.k) ')'], ...
                        'c^{C*}_{2,"S2"} (\sigma_2=0, k=0)', 'c_1', ...
                        'c^{C*}_{2,"S1"} (\sigma_2=0, k=0)', ...
                        ['c^{C*}_{2,"S1"} (\sigma_2=' num2str(t2model.sd2) ', k=' num2str(t2model.k) ')']);
            case 3
                title(['c_1 = ' num2str(c)]);
            case 4
                ylabel(['p(correct_1)_T = ' num2str(pct)]);
        end
    end
end