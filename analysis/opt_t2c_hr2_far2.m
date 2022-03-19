% opt_t2c_hr2_far2
%
% this script investigates the type 2 criteria that optimize HR2 - FAR2 
% according to the classical SDT model by plotting the effect of 
% d' and c on the optimal criteria.
%
% this script mirrors the analaysis of Figure 7 in the main manuscript.
%
% the details of this analysis can be readily changed by appropriately 
% editing the "parameter settings" section.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

clear

addpath(genpath('..'));

%% parameter settings

cs   = [-1, 0, 1];
ds   = .5 : .1 : 3;


%% construct the plot

figure;

i_plot = 0;
for i_c = 1:3

    c   = cs(i_c);

    % calculate the optimal type 2 criteria
    out = opt_t2c('hr2-far2', ds, c);

    i_plot = i_plot + 1;

    subplot(1, 3, i_plot); hold on;
    plot(ds, out.c2_rS2, 'r-', 'LineWidth', 1);
    plot(ds, cs(i_c)*ones(size(ds)), 'k-')
    plot(ds, out.c2_rS1, 'r--', 'LineWidth', 1);
    plot(ds, zeros(size(ds)), 'k--')

    xlabel('d''')
    xlim([min(ds), max(ds)])
    ylim([-3, 3])

    title(['c_1 = ' num2str(c)]);
    
    if i_plot == 2
        legend('c^{HF*}_{2,"S2"}', 'c_1', 'c^{HF*}_{2,"S1"}');
    end
end