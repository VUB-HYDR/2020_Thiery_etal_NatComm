

% --------------------------------------------------------------------
% function to plot LR as bar plots
% --------------------------------------------------------------------


function [] = mf_plot_LR(LR_bars, percentages, ylims, panel, region, flag_plot_legend)



% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% line colors
colors    = mf_colors;
colors_LR = [colors(32,:); ...
             colors(30,:); ...
             colors(39,:)     ];

                               
% define axes color                               
axcolor = [0.3 0.3 0.3]; % 70% contrast (so 0.3) is advised
         
         

% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% get number of probabilities(=percentages)
npercentages = length(percentages);


% get x-axis ticklabels
for i=1:npercentages
    xticklabels{i} = ['>' num2str(percentages(i)) '%']; %#ok<SAGROW>
end


% rescale for plotting
[LR_plot] = mf_riskratio_plotscale(LR_bars);



% --------------------------------------------------------------------
% Visualisation
% --------------------------------------------------------------------


% design figure
figure
set(gcf, 'color', 'w');


% plot data
h = bar(1:npercentages, LR_plot); hold on;
set(h(1),'BaseValue',1)
set(h,'edgecolor','none');
set(h(1),'facecolor', colors_LR(1,:));
set(h(2),'facecolor', colors_LR(2,:));
set(h(3),'facecolor', colors_LR(3,:));
axis([0.5 npercentages+0.5 mf_riskratio_plotscale(ylims(1,1)) mf_riskratio_plotscale(ylims(1,2))]);
set(gca, 'Fontsize', 15, 'Fontweight', 'Bold', 'Xcolor', axcolor, 'Ycolor', axcolor);
set(gca, 'xtick', 1:npercentages, 'xticklabel', xticklabels, 'XTickLabelRotation', 90);
set(gca, 'ytick', mf_riskratio_plotscale([1/32 1/16 1/8 1/6 1/5 1/4 1/3 1/2 1 2 3 4 5 6 8 16 32]), 'YTickLabel', {'1/32','1/16','1/8','1/6','1/5','1/4','1/3','1/2','1','2','3','4','5','6','8','16','32'})
% xlabel('T_{2m,max} percentiles', 'Fontsize', 18, 'Fontweight', 'Bold'); 
ylabel('LR', 'Fontsize', 18, 'Fontweight', 'Bold'); 


% plot legend
if flag_plot_legend == 1
    l(1) = plot(NaN,1, 'color', colors_LR(1,:), 'linewidth', 3); % dummy plot for legend
    l(2) = plot(NaN,1, 'color', colors_LR(2,:), 'linewidth', 3); % dummy plot for legend
    l(3) = plot(NaN,1, 'color', colors_LR(3,:), 'linewidth', 3); % dummy plot for legend
    % legend(l, 'present-day warming','present-day irrigation','present-day warming + irrigation', 'location', 'southwest');
    legend(l, 'All forcings except irrigation','Irrigation expansion','All forcings', 'location', 'southwest');
    set(legend, 'box', 'off', 'Fontweight', 'Bold', 'Fontsize', 14, 'textcolor', axcolor);
end

% add text
text(0.6,ylims(1,2),panel,'ver','bottom','hor','center','Fontsize', 18, 'color', axcolor)
text(npercentages+0.5,ylims(1,2), region,'ver','bottom','hor','right','Fontsize', 18, 'Fontweight', 'Bold', 'color', axcolor)


end
