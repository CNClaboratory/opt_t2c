function [opt_sim, opt_calc, M] = opt_t2c_sim(type, param, t2model, settings)
% [opt_sim, opt_calc, M] = opt_t2c_sim(type, param, t2model, settings)
%
% Perform SDT simulations to find the optimal type 2 criteria for a given
% optimization objective, given the type 1 and type 2 parameters of the
% model.
%
% INPUTS
% ------
% * type - a string that determines what target measure is to be optimized. 
%     Possible values of type are:
%     * 't2acc'       - p(type 2 correct)
%     * 't2reward'    - expected reward based on type 2 reward contingency table
%     * 'calibration' - reporting high confidence for trials where p(correct) >
%                       some p(correct) threshold
%     * 'hr2-far2'    - difference b/t type 2 hit rate and type 2 false alarm rate
%
% * param - a struct holding the main model parameters.
% 
%     Required fields of param are:
%     param.d   - type 1 sensitivity d'
%     param.c   - type 1 criterion c
%     param.pS2 - probability of S2 stimulus presentation
%
%     When type is 't2reward', param.R must be specified as follows:
%     param.R   - a struct with fields R.hit2, R.miss2, R.CR2, R.FA2, where each
%         field contains the amount of points earned following each type 2 event
%         (negative values indicate loss of points and zero values indicate no
%         points won or lost).
%
%     When type is 'calibration', param.pcorr_thresh must be specified as follows:
%     param.pcorr_thresh - p(correct)_T, i.e. the threshold level of p(correct) 
%         such that confidence is high when estimated p(correct) > p(correct)_T
%
%     if R or pcorr_thresh are included as fields of param for
%     non-applicable types, these fields are ignored.
%
% * t2model - a struct holding details of the type 2 model.
%
%     Required field is:
%     t2model.t2fn - a function handle for the type 2 function used to
%         transform a 1 x ntrials input vector of type 1 evidence samples x1 
%         into a 1 x ntrials output vector of type 2 evidence samples x2.
%         This function should be formatted as
%         
%         x2 = t2fn(x1, param);
%
%         where param is a struct containing fields for the parameters used
%         by t2fn to compute x2. 
%
%         For example, see the function type2noisySignalLoss included in this 
%         toolbox. If using this function, the syntax would be
% 
%         t2model.t2fn = @type2noisySignalLoss;
%
%     t2model should also contain fields corresponding to the parameters of
%     the model used by t2fn. For example, the included function type2noisySignalLoss 
%     takes an input param with fields k and sd2 specifying the signal loss
%     parameter and the type 2 noise parameter, respectively. Thus the
%     t2model struct must also contain fields for k and sd2 specifying the
%     values to be used in the simulation.
%
%     t2model can also contain an optional field t2titlestr, which is used
%     to provide details of the type 2 model in the final line of the title of 
%     the plot that is optionally produced by opt_t2c_sim.
%
% * settings - a struct holding details of the simulation settings.
%
%     Required fields of settings are:
%     settings.ntrials       - the number of trials to simulate
%     settings.compute_metad - set to 1 if you wish the simulation results
%         to include meta-d' analysis, 0 otherwise. setting to 1 can
%         substantially increase computation time.
%     settings.makePlot      - set to 1 if you wish to make a plot of the 
%         simulation results, 0 otherwise.
%
% OUTPUTS
% -------
% * opt_sim - a struct holding the results of the simulation.
%  
%     Fields of the struct are:
%     * c2_rS1     - the simulation-derived optimal type 2 criterion for "S1" responses.
%     * c2_rS2     - the simulation-derived optimal type 2 criterion for "S2" responses.
%     * t2perf_rS1 - the simulated type 2 performance measure at c2_rS1.
%     * t2perf_rS2 - the simulated type 2 performance measure at c2_rS2.
%
%     "Type 2 performance" refers to the measure to be optimized. This is
%     straightforward to interpret for the cases of maximizing type 2 accuracy, 
%     type 2 reward, or HR2 - FAR2. For instance, for type 2 reward,
%     t2perf_rS1 denotes the average reward accrued for "S1" responses, and
%     so on.
%
%     For calibration, "type 2 performance" refers to the ideal observer's estimated 
%     p(correct) at the type 2 criterion, where the type 2 criterion is
%     selected so as to make estimated p(correct) as close as possible to
%     the target value of p(correct)_T. (See XXX for further discussion of
%     how calibration is simulated.)
% 
% * opt_calc - a struct holding the optimal type 2 criteria according to
%     equations derived under classical SDT (which assumes that type 2
%     evidence samples are identical to type 1 evidence samples).
%
%     Fields of the struct include:
%     * c2_rS1     - the theory-derived optimal type 2 criterion for "S1"
%                    responses, according to classical SDT.
%     * c2_rS2     - the theory-derived optimal type 2 criterion for "S2"
%                    responses, according to classical SDT.
%     * t2perf_rS1 - the *simulated* type 2 performance measure evaluated at 
%                    the *simulated* c2_rS1 closest to the theoretically
%                    derived c2_rS1 value.
%     * t2perf_rS2 - the *simulated* type 2 performance measure evaluated at 
%                    the *simulated* c2_rS2 closest to the theoretically
%                    derived c2_rS2 value.
%
%     Remaining fields in this struct are as described in the documentation
%     for the opt_t2c function.
%
% * M - a struct containing results of meta-d' analysis. This analysis is
%     conducted by setting 9 evenly spaced type 2 criteria on either side
%     of the type 1 criterion, comparing type 2 evidence samples to these
%     criteria to derive confidence ratings for each simulated trial, and
%     using these confidence ratings to conduct meta-d' analysis on the
%     simulated data.
%
%     Fields of the struct are:
%     * md     - overall meta-d'
%     * M      - overall M_ratio
%     * md_rS1 - meta-d' for "S1" responses
%     * md_rS2 - meta-d' for "S2" responses
%     * M_rS1  - M_ratio for "S1" responses
%     * M_rS2  - M_ratio for "S2" responses
%
%     if settings.compute_metad == 0, then the fields of M are set to NaN.
%
% 3/19/2022  Brian Maniscalco, Lucie Charles, & Megan Peters


%% parse inputs & prepare for analysis

% add toolbox folder to path
addpath(genpath('toolbox'));

% unpack inputs
v2struct(param);
v2struct(t2model);
v2struct(settings);

% define search space for type 2 criteria
c2_rS1 = c : -.01 : c-4;
c2_rS2 = c :  .01 : c+4;

pS1 = 1 - pS2;

%% calculate optimal type 2 criteria under SDT assumptions 
%  (i.e. assuming x2 == x1)

switch type
    
    % type 2 accuracy
    case 't2acc'
        opt_calc  = opt_t2c(type, d, c, pS2);
        legendstr = 'p(correct_2)';

    % type 2 reward
    case 't2reward'
        Q2 = (R.CR2 - R.FA2) / (R.hit2 - R.miss2);
        opt_calc = opt_t2c(type, d, c, pS2, Q2);
        legendstr = 'E(reward_2)';
        
    % calibrating confidence to expected accuracy
    case 'calibration'
        OT = pcorr_thresh / (1-pcorr_thresh);
        opt_calc  = opt_t2c(type, d, c, pS2, OT);
        legendstr = 'p(correct_1)';
        
    % type 2 HR - FAR
    case 'hr2-far2'
        opt_calc = opt_t2c(type, d, c);
        legendstr = 'HR_2 - FAR_2';

end


%% conduct the simulation for type 1 performance

% compute # of trials for each stim type
nS1 = round(pS1 * ntrials);
nS2 = ntrials - nS1;

% simulate type 1 and type 2 evidence values
x1 = [ normrnd(-d/2, 1, [1,nS1]), normrnd(+d/2, 1, [1,nS2]) ]; % evidence samples for type 1 decisions
x2 = t2fn(x1, t2model);                                        % evidence samples for type 2 decisions

% compute stimID, response, and accuracy
stimID = [zeros(1, nS1), ones(1, nS2)];  % 0 -> stim=S1,   1 -> stim=S2 
resp   = x1 > c;                         % 0 -> resp="S1", 1 -> resp="S2"
acc    = resp == stimID;                 % 0 -> incorrect, 1 -> correct


%% conduct the simulation for type 2 performance

% initialize arrays for speed
t2perf_rS1 = nan(size(c2_rS1));
t2perf_rS2 = nan(size(c2_rS1));

% compute type 2 performance for all values in the type 2 criteria search space 
for i = 1:length(c2_rS2)

    % get confidence on each trial for current type 2 criteria
    conf_rS1 = (x1 <= c) & (x2 < c2_rS1(i)); % 1 for high conf "S1" trials, 0 otherwise
    conf_rS2 = (x1 >  c) & (x2 > c2_rS2(i)); % 1 for high conf "S2" trials, 0 otherwise
    
    % compute the measure of type 2 performance
    switch type

        % type 2 accuracy
        case 't2acc'
            % p(acc==conf | resp="S1")
            t2perf_rS1(i) = ( sum(resp==0 & acc==1 & conf_rS1==1) + ... % total # of type 2 hits
                              sum(resp==0 & acc==0 & conf_rS1==0) ) ... % total # of type 2 CRs
                              / sum(resp==0);                           % total # of "S1" responses
            
            % p(acc==conf | resp="S2")
            t2perf_rS2(i) = ( sum(resp==1 & acc==1 & conf_rS2==1) + ... % total # of type 2 hits
                              sum(resp==1 & acc==0 & conf_rS2==0) ) ... % total # of type 2 CRs
                              / sum(resp==1);                           % total # of "S2" responses
            
            
        % type 2 reward
        case 't2reward'
            % points per trial | resp = "S1"
            t2perf_rS1(i) = ( R.hit2  * sum(resp==0 & acc==1 & conf_rS1==1) + ... % points from type 2 hits
                              R.miss2 * sum(resp==0 & acc==1 & conf_rS1==0) + ... % points from type 2 misses
                              R.FA2   * sum(resp==0 & acc==0 & conf_rS1==1) + ... % points from type 2 FAs
                              R.CR2   * sum(resp==0 & acc==0 & conf_rS1==0) ) ... % points from type 2 CRs
                              / sum(resp==0);                                     % total # of "S1" responses
                          
            % points per trial | resp = "S2"
            t2perf_rS2(i) = ( R.hit2  * sum(resp==1 & acc==1 & conf_rS2==1) + ... % points from type 2 hits
                              R.miss2 * sum(resp==1 & acc==1 & conf_rS2==0) + ... % points from type 2 misses
                              R.FA2   * sum(resp==1 & acc==0 & conf_rS2==1) + ... % points from type 2 FAs
                              R.CR2   * sum(resp==1 & acc==0 & conf_rS2==0) ) ... % points from type 2 CRs
                              / sum(resp==1);                                     % total # of "S2" responses
                      
                          
        % calibrating confidence to expected accuracy
        case 'calibration'
            bin_width = .01;
            
            % define the evidence type used for calibration.
            % in some sense it is most intuitive for the observer to use 
            % the type 1 evidence variable x1 for calibration, since 
            % calibration requires evaluating p(correct) and p(correct) 
            % depends on x1. however, this implies that calibration is 
            % independent of type 2 evidence x2, which contrasts with the 
            % other types of optimization contexts. thus here we default to 
            % setting x2 as the type of evidence used for calibration, in 
            % spite of this being a somewhat counterintuitive choice. 
            x = x2;
            
            % compute estimated p(correct) when evidence samples are 
            % contained within a small bin centered on c2_rS1. the estimate
            % for "S1" responses is obtained by finding the proportion of 
            % evidence samples in this bin that originate from S1 stimuli
            x_rS1_filter  = (x >  c2_rS1(i) - bin_width/2) & ...
                            (x <= c2_rS1(i) + bin_width/2);
            t2perf_rS1(i) = sum(x_rS1_filter & stimID==0) / sum(x_rS1_filter);
            
            % compute estimated p(correct) when evidence samples are 
            % contained within a small bin centered on c2_rS2. the estimate
            % for "S2" responses is obtained by finding the proportion of 
            % evidence samples in this bin that originate from S2 stimuli
            x_rS2_filter  = (x >  c2_rS2(i) - bin_width/2) & ...
                            (x <= c2_rS2(i) + bin_width/2);
            t2perf_rS2(i) = sum(x_rS2_filter & stimID==1) / sum(x_rS2_filter);

            
        % type 2 HR - FAR
        case 'hr2-far2'
            
            % compute type 2 HR and FAR
            HR2_rS1  = sum(resp==0 & acc==1 & conf_rS1==1) / sum(resp==0 & acc==1);
            FAR2_rS1 = sum(resp==0 & acc==0 & conf_rS1==1) / sum(resp==0 & acc==0);
            HR2_rS2  = sum(resp==1 & acc==1 & conf_rS2==1) / sum(resp==1 & acc==1);
            FAR2_rS2 = sum(resp==1 & acc==0 & conf_rS2==1) / sum(resp==1 & acc==0);
            
            % compute HR2-FAR2 for each response type
            t2perf_rS1(i) = HR2_rS1 - FAR2_rS1;
            t2perf_rS2(i) = HR2_rS2 - FAR2_rS2;
            
    end
end


%% optionally compute meta-d' using evenly distributed conf criteria

% default meta-d' output to nans
M.md     = nan;
M.M      = nan;
M.md_rS1 = nan;
M.md_rS2 = nan;
M.M_rS1  = nan;
M.M_rS2  = nan;

% define evenly distributed type 2 criteria
nRatings  = 10;                              % use 10-point rating scale to ensure robust estimation of meta-d'
c2_rS1_md = c - linspace(.1, 4, nRatings-1); % evenly distribute type 2 criteria for "S1" responses across a wide range of the decision axis
c2_rS2_md = c + linspace(.1, 4, nRatings-1); % evenly distribute type 2 criteria for "S2" responses across a wide range of the decision axis

if compute_metad

    % get confidence rating on each trial
    rating = ones(size(x1));
    for i = 1:length(c2_rS1_md)
        rating( x1 <= c & x2 < c2_rS1_md(i) ) = i + 1; % S1 responses
        rating( x1 >  c & x2 > c2_rS2_md(i) ) = i + 1; % S2 responses
    end
    
    % convert trial data to response counts
    [nR_S1, nR_S2] = trials2counts(stimID, resp, rating, nRatings, 1);

    % compute overall meta-d'
    fit  = fit_meta_d_MLE(nR_S1, nR_S2);
    M.md = fit.meta_da;
    M.M  = fit.M_ratio;

    % compute response-conditional meta-d'
    fit_rs   = fit_rs_meta_d_MLE(nR_S1, nR_S2);
    M.md_rS1 = fit_rs.meta_da_rS1;
    M.md_rS2 = fit_rs.meta_da_rS2;
    M.M_rS1  = fit_rs.M_ratio_rS1;
    M.M_rS2  = fit_rs.M_ratio_rS2;
end


%% find best-performing type 2 criteria in the simulation

switch type

    % for calibration, optimal type 2 criteria minimize the magnitude of 
    % the difference between pcorr and pcorr_thresh
    case 'calibration'
        pc_rS1_diff        = abs( t2perf_rS1 - pcorr_thresh );
        [m, ind]           = min(pc_rS1_diff);
        opt_sim.c2_rS1     = c2_rS1(ind);
        opt_sim.t2perf_rS1 = t2perf_rS1(ind);
        
        pc_rS2_diff        = abs( t2perf_rS2 - pcorr_thresh );
        [m, ind]           = min(pc_rS2_diff);
        opt_sim.c2_rS2     = c2_rS2(ind);
        opt_sim.t2perf_rS2 = t2perf_rS2(ind);
        
        
    % for t2acc, t2reward, and hr2-far2, optimal type 2 criteria maximize
    % the type 2 performance measure
    otherwise
        [m, ind]           = max(t2perf_rS1);
        opt_sim.c2_rS1     = c2_rS1(ind);
        opt_sim.t2perf_rS1 = t2perf_rS1(ind);

        [m, ind]           = max(t2perf_rS2);
        opt_sim.c2_rS2     = c2_rS2(ind);
        opt_sim.t2perf_rS2 = t2perf_rS2(ind);
        
end


%% find simulated type 2 performance corresponding to type 2 criteria calculated under SDT assumptions
%  i.e. if the observer assumes x2==x1 and sets their type 2 criteria
%  optimally according to that assumption, what is their resulting type 2
%  performance?

c2_rS1_diff         = abs( c2_rS1 - opt_calc.c2_rS1 );
[m, ind]            = min(c2_rS1_diff);
opt_calc.t2perf_rS1 = t2perf_rS1(ind);

c2_rS2_diff         = abs( c2_rS2 - opt_calc.c2_rS2 );
[m, ind]            = min(c2_rS2_diff);
opt_calc.t2perf_rS2 = t2perf_rS2(ind);


%% plot simulation results
% plots type 2 performance as a function of type 2 criterion for each
% response type and displays the optimal type 2 criteria derived from
% simulation results alongside the theoretical type 2 criteria which ignore
% sources of type 2 noise

if makePlot
    
    c_minmax = [0, max([t2perf_rS1,t2perf_rS2])];
    
    figure;
    
    % plot "S1" response results
    subplot(1,2,1); hold on;
    plot(c*[1,1], c_minmax, 'k-', 'LineWidth', 2)
    plot(c2_rS1, t2perf_rS1); 
    plot(opt_sim.c2_rS1*[1,1],  c_minmax, 'r--', 'LineWidth', 2)
    plot(opt_calc.c2_rS1*[1,1], c_minmax, 'b:', 'LineWidth', 2)

    if strcmp(type, 'calibration')
        plot(c2_rS1, pcorr_thresh*ones(size(c2_rS1)), 'g-')
        legend('c_1', [legendstr ' at c_2 (sim)'], 'opt c_2 accounting for type 2 model (sim)', 'opt c_2 ignoring type 2 model (calc)', 'p(correct)_T', 'location', 'southwest')    
    else
        legend('c_1', [legendstr ' at c_2 (sim)'], 'opt c_2 accounting for type 2 model (sim)', 'opt c_2 ignoring type 2 model (calc)', 'location', 'southwest')    
    end
    xlabel('c_2 for "S1" responses')

    
    % plot "S2" response results
    subplot(1,2,2); hold on;
    plot(c*[1,1], c_minmax, 'k-', 'LineWidth', 2)
    plot(c2_rS2, t2perf_rS2); 
    plot(opt_sim.c2_rS2*[1,1],  c_minmax, 'r--', 'LineWidth', 2)
    plot(opt_calc.c2_rS2*[1,1], c_minmax, 'b:', 'LineWidth', 2)

    if strcmp(type, 'calibration')
        plot(c2_rS2, pcorr_thresh*ones(size(c2_rS2)), 'g-')
    end
    xlabel('c_2 for "S2" responses')

    
    % title
    switch type
        case 't2reward'
            titlestr_extra = [', Q_2 = ' num2str(Q2)];
        case 'calibration'
            titlestr_extra = [', O_T = ' num2str(OT)];
        otherwise
            titlestr_extra = [];
    end  

    titlestr = {['simulation results for ' type ' (n trials = ' num2str(ntrials) ')'], ...
                ['with d'' = ' num2str(d) ', c_1 = ' num2str(c) ', p(S2) = ' num2str(pS2), titlestr_extra], ...
                 t2titlestr};

    try
        sgtitle(titlestr);
    catch
        title(titlestr);
    end
end

end