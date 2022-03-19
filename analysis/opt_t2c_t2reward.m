% opt_t2c_t2reward
%
% this script investigates the type 2 criteria that optimize type 2
% reward according to the classical SDT model by plotting the effect of 
% d', c, and Q2 on the optimal criteria.
%
% this script mirrors the analaysis of Figure 5 in the main manuscript.
%
% the details of this analysis can be readily changed by appropriately 
% editing the "parameter settings" section.
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

clear

addpath(genpath('..'));

%% parameter settings

pS2 = 0.5;
Q2s = [1/3, 1, 3];
cs  = [-1, 0, 1];
ds  = .5 : .1 : 3;


%% construct the plot

figure;

i_plot = 0;
for i_Q2 = 1:3
    for i_c = 1:3
        
        c  = cs(i_c);
        Q2 = Q2s(i_Q2);
        
        % calculate the optimal type 2 criteria
        out = opt_t2c('t2reward', ds, c, pS2, Q2);
        
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
                ylabel(['Q_2 = ' num2str(Q2)]);
            case 2
                title(['c_1 = ' num2str(c)]);
                legend('c^{R*}_{2,"S2"}', 'c_1', 'c^{R*}_{2,"S1"}');
            case 3
                title(['c_1 = ' num2str(c)]);
            case 4
                ylabel(['Q_2 = ' num2str(Q2)]);
            case 7
                ylabel(['Q_2 = ' num2str(Q2)]);                
        end
    end
end