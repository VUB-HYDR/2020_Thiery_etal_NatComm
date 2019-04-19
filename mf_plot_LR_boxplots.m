

% --------------------------------------------------------------------
% function to plot LR as box plots
% --------------------------------------------------------------------


function [] = mf_plot_LR_boxplots(LR_bars, boxplotdata_LR_pic2ctl, boxplotdata_LR_ctl2irr, boxplotdata_LR_pic2irr, percentages, ylims, panel, region, flag_plot_legend)



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


% initialise boxplot offset
boxplot_offset = 0.15;
       


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


% construct matrices for plotting and rescale
[boxplotdata_LR_pic2ctl_pinb_plot] = mf_riskratio_plotscale(boxplotdata_LR_pic2ctl);
[boxplotdata_LR_ctl2irr_pinb_plot] = mf_riskratio_plotscale(boxplotdata_LR_ctl2irr);
[boxplotdata_LR_pic2irr_pinb_plot] = mf_riskratio_plotscale(boxplotdata_LR_pic2irr);



% --------------------------------------------------------------------
% Visualisation
% --------------------------------------------------------------------


% design figure
figure
set(gcf, 'color', 'w');


% plot dummy data
h = bar(1:npercentages, LR_plot .*0 + 1); hold on;
set(h(1),'BaseValue',1)
set(h,'edgecolor','none');
colormap(colors_LR)

axis([0.5 npercentages+0.5 mf_riskratio_plotscale(ylims(1,1)) mf_riskratio_plotscale(ylims(1,2))]);
set(gca,'xtickmode','auto','xticklabelmode','auto');                               % remove ticklabels so that originals appear

set(gca, 'Fontsize', 15, 'Fontweight', 'Bold', 'Xcolor', axcolor, 'Ycolor', axcolor, 'xtick', 1:npercentages, 'xticklabel', xticklabels);
text_h = findobj(gca, 'Type', 'text');
for i = 1:length(text_h)
    set(text_h(i), 'FontSize', 11, 'Fontweight', 'Bold', 'color', axcolor, 'Rotation', 90, 'String', xticklabels{length(xticklabels)-i+1}, 'HorizontalAlignment', 'right')
end

% smaller box for axes, in order to un-hide the labels
squeeze = 0.05;
left    = 0.02;
right   = 1;
bottom  = squeeze;
top     = 1-squeeze;
set(gca, 'OuterPosition', [left bottom right top])

% plot actual boxplots
boxplot(boxplotdata_LR_pic2ctl_pinb_plot, 'PlotStyle', 'compact', 'colors', colors_LR(1,:), 'symbol','w.', 'positions', (1:npercentages)-boxplot_offset , 'labels', repmat({''},1,npercentages)); hold on; % plot boxplots
boxplot(boxplotdata_LR_ctl2irr_pinb_plot, 'PlotStyle', 'compact', 'colors', colors_LR(2,:), 'symbol','w.', 'positions', (1:npercentages)                , 'labels', repmat({''},1,npercentages)); hold on; % plot boxplots
boxplot(boxplotdata_LR_pic2irr_pinb_plot, 'PlotStyle', 'compact', 'colors', colors_LR(3,:), 'symbol','w.', 'positions', (1:npercentages)+boxplot_offset , 'labels', repmat({''},1,npercentages)); hold on; % plot boxplots
axis([0.5 npercentages+0.5 mf_riskratio_plotscale(ylims(1,1)) mf_riskratio_plotscale(ylims(1,2))]);

set(gca, 'xtick', 1:npercentages, 'xticklabel', xticklabels, 'XTickLabelRotation', 90);
set(gca, 'Ytick', mf_riskratio_plotscale([1/32 1/16 1/8 1/6 1/5 1/4 1/3 1/2 1 2 3 4 5 6 8 16 32]),'YTickLabel', {'1/32','1/16','1/8','1/6','1/5','1/4','1/3','1/2','1','2','3','4','5','6','8','16','32'})
% set(gca,'YTickLabel',{'1/64','1/32','1/16','1/8','1','8','16'}, 'Ytick', mf_riskratio_plotscale([1/64 1/32 1/16 1/8 1 8 16]))
% xlabel('TX percentiles', 'Fontsize', 18, 'Fontweight', 'Bold'); 
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
