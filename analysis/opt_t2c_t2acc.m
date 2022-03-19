% opt_t2c_t2acc
%
% this script investigates the type 2 criteria that optimize type 2
% accuracy according to the classical SDT model by plotting the effect of 
% d', c, and p(S2) on the optimal criteria.
%
% this script mirrors the analaysis of Figure 4 in the main manuscript.
%
% the details of this analysis can be readily changed by appropriately 
% editing the "parameter settings" section.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

clear

addpath(genpath('..'));

%% parameter settings

pS2s = [.5, .3, .15];
cs   = [-1, 0, 1];
ds   = .5 : .1 : 3;


%% construct the plot

figure;

i_plot = 0;
for i_pS2 = 1:3
    for i_c = 1:3
        
        c   = cs(i_c);
        pS2 = pS2s(i_pS2);
        
        % calculate the optimal type 2 criteria
        out = opt_t2c('t2acc', ds, c, pS2);
        
        i_plot = i_plot + 1;
        
        subplot(3, 3, i_plot); hold on;
        plot(ds, out.c2_rS2, 'r-', 'LineWidth', 1);
        plot(ds, cs(i_c)*ones(size(ds)), 'k-')
        plot(ds, out.c2_rS1, 'r--', 'LineWidth', 1);
        plot(ds, zeros(size(ds)), 'k--')
        
        xlabel('d''')
        xlim([min(ds), max(ds)])
        ylim([-3, 3])
        
        switch i_plot
            case 1
                title(['c_1 = ' num2str(c)]);
                ylabel(['p(S2) = ' num2str(pS2)]);
            case 2
                title(['c_1 = ' num2str(c)]);
                legend('c^{A*}_{2,"S2"}', 'c_1', 'c^{A*}_{2,"S1"}');
            case 3
                title(['c_1 = ' num2str(c)]);
            case 4
                ylabel(['p(S2) = ' num2str(pS2)]);
            case 7
                ylabel(['p(S2) = ' num2str(pS2)]);                
        end
    end
end