
% --------------------------------------------------------------------
% visualisation subroutine
% note: preferably run "main"
% --------------------------------------------------------------------


% clean up
clc;
close all;


% flags for paper plots
flags.plot_fig1   = 0; % 0: do not plot figure 1 of paper
                       % 1: plot figure 1 of paper
flags.plot_fig2   = 0; % 0: do not plot figure 2 of paper
                       % 1: plot figure 2 of paper
flags.plot_fig3   = 0; % 0: do not plot figure 3 of paper
                       % 1: plot figure 3 of paper
flags.plot_sfig1  = 0; % 0: do not plot supplementary figure 1 of paper
                       % 1: plot supplementary figure 1 of paper
flags.plot_sfig2  = 0; % 0: do not plot supplementary figure 2 of paper
                       % 1: plot supplementary figure 2 of paper
flags.plot_sfig3  = 0; % 0: do not plot supplementary figure 3 of paper
                       % 1: plot supplementary figure 3 of paper
flags.plot_sfig4  = 1; % 0: do not plot supplementary figure 4 of paper
                       % 1: plot supplementary figure 4 of paper
flags.plot_sfig5  = 0; % 0: do not plot supplementary figure 5 of paper
                       % 1: plot supplementary figure 5 of paper
flags.plot_sfig6  = 0; % 0: do not plot supplementary figure 6 of paper
                       % 1: plot supplementary figure 6 of paper
flags.plot_sfig7  = 0; % 0: do not plot supplementary figure 7 of paper
                       % 1: plot supplementary figure 7 of paper
flags.plot_sfig8  = 0; % 0: do not plot supplementary figure 8 of paper
                       % 1: plot supplementary figure 8 of paper
flags.plot_sfig9  = 0; % 0: do not plot supplementary figure 9 of paper
                       % 1: plot supplementary figure 9 of paper

                     
% flags for other plots               
flags.plot_surv  = 0; % 0: do not plot survival plots
                      % 1: plot survival plots
flags.plot_LRts  = 0; % 0: do not plot Risk Ratios (ts: LR computed on area-averaged time series)
                      % 1: plot Risk Ratios (ts: LR computed on area-averaged time series)
        

                     
% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% set colorscale axes
caxes.LR   = [1/8         8  ];
caxes.Pirr = [  0        50  ];
caxes.popd = [  0       500  ];
caxes.regr = [ -1.5       1.5];
caxes.mon  = [  0.5        12.5  ];
caxes.dmon = [ -3.5       3.5];


% set colormaps
colormaps.LR   =        mf_colormap_cpt('dkbluered'        , 28);
colormaps.Pirr =        mf_colormap_cpt('cbacYlGnBu09'     , 10);
colormaps.popd =        mf_colormap_cpt('RdPu_09'          , 9 );
colormaps.regr =        mf_colormap_cpt('dkbluered'        , 30);
colormaps.mon  =        mf_colormap_cpt('passionfruit'     , 12);
% colormaps.mon  =        mf_colormap_cpt('Paired_12'     , 12);
colormaps.dmon =        mf_colormap_cpt('cbcBrBG'          ,  7);


% Set center of change to white
colormaps.LR(14:15,:)   = [1 1 1; 1 1 1];
colormaps.regr(15:16,:) = [1 1 1; 1 1 1];
colormaps.popd(10,:)    = [0.15  0  0.30];
colormaps.popd(1,:)     = [1     1  1   ];

colormaps.Pirr(1:9,:)   = colormaps.Pirr(2:10,:);
colormaps.Pirr(10,:)    = [0.01 0.07 0.23];

colormaps.mon(1,:)   = [0.96 0.96 0];
colormaps.mon(end,:) = [0 0 0];

% line colors
colors = mf_colors;

                               
% define axes color                               
axcolor = [0.3 0.3 0.3]; % 70% contrast (so 0.3) is advised


% get alphabet                               
alphabet = char(repmat('a' + (1:26) - 1, [2 1]))';



% --------------------------------------------------------------------
% paper figure 1
% --------------------------------------------------------------------


if flags.plot_fig1 == 1

    
disp('see python scripts from Auke Visser, reworked by Wim Thiery and Matthias Hauser');
disp('/net/exo/landclim/mathause/analysis/2018_for_wim_irrig/scripts/CRU_CESM_binnedboxplot_regional_WT_MH.ipynb');
disp('Mathias Huaser: I did it in a jupyter notebook - by executing ./startup it should launch the correct environment.');


end



% --------------------------------------------------------------------
% paper figure 2
% --------------------------------------------------------------------


if flags.plot_fig2 == 1 && flags.refr == 0


% rescale for plotting
LR_pic2ctl_plot = mf_riskratio_plotscale(LR_pic2ctl);
LR_ctl2irr_plot = mf_riskratio_plotscale(LR_ctl2irr);
LR_pic2irr_plot = mf_riskratio_plotscale(LR_pic2irr);
caxes.LR_plot   = mf_riskratio_plotscale(caxes.LR);


% loop over percentages and plot LR maps
% for i=1:npercentages
for i=5  % P99 only

    % pic2ctl LR map
    cbh = mf_plot_dom2(lon_mod, lat_mod, LR_pic2ctl_plot(:,:,i), [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'a', 'All forcings except irrigation', 'LR'); hold on;
    % island = island & ispinb;
    % cbh = mf_plot_dom2(lon_mod, lat_mod, LR_pic2ctl_plot(:,:,i) .* ispinb, [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'a', ['P' num2str(percentages(i)) ' - pic2ctl'], 'LR'); hold on;
    set(cbh,'XTickLabel',{'1/8','1/6','1/4','1/2','1','2','4','6','8'}, 'Xtick', mf_riskratio_plotscale([1/8 1/6 1/4 1/2 1 2 4 6 8]))
    export_fig(sprintf('figures/LR_map_pic2ctl_P%s_lowres.png' , num2str(percentages(i),'%.3f')), '-transparent');
%     export_fig(sprintf('figures/LR_map_pic2ctl_P%s_highres.png' , num2str(percentages(i),'%.3f')), '-transparent', '-m10'); close all;


    % ctl2irr LR map
    cbh = mf_plot_dom2(lon_mod, lat_mod, LR_ctl2irr_plot(:,:,i), [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'b', 'Irrigation expansion', 'LR'); hold on;
    % island = island .* 0;
    % cbh = mf_plot_dom2(lon_mod, lat_mod, LR_ctl2irr_plot(:,:,i) .* ispinb, [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'b', ['P' num2str(percentages(i)) ' - ctl2irr'], 'LR'); hold on;
    set(cbh,'XTickLabel',{'1/8','1/6','1/4','1/2','1','2','4','6','8'}, 'Xtick', mf_riskratio_plotscale([1/8 1/6 1/4 1/2 1 2 4 6 8]))
    export_fig(sprintf('figures/LR_map_ctl2irr_P%s_lowres.png' , num2str(percentages(i),'%.3f')), '-transparent');
%     export_fig(sprintf('figures/LR_map_ctl2irr_P%s_highres.png', num2str(percentages(i),'%.3f')), '-transparent', '-m10'); close all;


    % pic2irr LR map
    cbh = mf_plot_dom2(lon_mod, lat_mod, LR_pic2irr_plot(:,:,i), [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'c', 'All forcings', 'LR'); hold on;
    % island = island .* 0;
    % cbh    = mf_plot_dom2(lon_mod, lat_mod, LR_pic2irr_plot(:,:,i) .* ispinb, [], island, caxes.LR_plot, colormaps.LR, 0, 2, 'c', ['P' num2str(percentages(i)) ' - pic2irr'], 'LR'); hold on;
    set(cbh,'XTickLabel',{'1/8','1/6','1/4','1/2','1','2','4','6','8'}, 'Xtick', mf_riskratio_plotscale([1/8 1/6 1/4 1/2 1 2 4 6 8]))
    export_fig(sprintf('figures/LR_map_pic2irr_P%s_lowres.png' , num2str(percentages(i),'%.3f')), '-transparent');
%     export_fig(sprintf('figures/LR_map_pic2irr_P%s_highres.png', num2str(percentages(i),'%.3f')), '-transparent', '-m10'); close all;
    
    
end
    

end



% --------------------------------------------------------------------
% paper figure 3
% --------------------------------------------------------------------


if flags.plot_fig3 == 1 && flags.refr == 0


% prepare LR bar plots
LR_bars_lp   = [LR_pic2ctl_lp LR_ctl2irr_lp LR_pic2irr_lp];
LR_bars_ip   = [LR_pic2ctl_ip LR_ctl2irr_ip LR_pic2irr_ip];
LR_bars_pinb = [LR_pic2ctl_pinb LR_ctl2irr_pinb LR_pic2irr_pinb];


% plot all  land
mf_plot_LR(LR_bars_lp, percentages, [1/6.25 6.25], 'a', 'All land', 1)
export_fig('figures/LR_median_lp.png', '-transparent'); % save figure


% plot irrigated land
mf_plot_LR(LR_bars_ip, percentages, [1/6.25 6.25], 'b', 'Irrigated land', 0)
export_fig('figures/LR_median_ip.png', '-transparent'); % save figure


% % % uncomment these for plotting only one bar (for conference presentations)
% % LR_bars_pinb(:,2:3) = 1;
% % LR_bars_pinb([1:4 6:end],1) = 1;


% plot PINB
mf_plot_LR(LR_bars_pinb, percentages, [1/6.25 6.25], 'c', 'South Asia', 0)
export_fig('figures/LR_median_pinb.png', '-transparent'); % save figure


end



% --------------------------------------------------------------------
% paper supplementary figure 1 - Area Equipped for Irrigation (AEI)
% --------------------------------------------------------------------


if flags.plot_sfig1 == 1


% define masks
isunirr_1915 = island & ~(pct_irr_1915 < 1);
isunirr      = island & ~(pct_irr < 1);
isinhabited  = island & ~(popdens < 1);


% percentage of pixel equipped for irrigation map - year 1915
cbh = mf_plot_dom2(lon_mod, lat_mod, pct_irr_1915, [], isunirr_1915, caxes.Pirr, colormaps.Pirr, 0, 2, 'a', 'Area equipped for irrigation - 1915', 'Grid cell fraction equipped for irrigation [%]'); hold on;
set(cbh,'YTick',[0:10:50],'YTicklabel',[1 10:10:50])
for i=1:length(pol_lat) % plot srex polygons and their names (repeat first corner at the end to close box)
    m_line([pol_lon{i}; pol_lon{i}(1)], [pol_lat{i}; pol_lat{i}(1)],'linewi',1,'color','r');    % Area outline
end
m_text(pol_lon{1}(2)   ,pol_lat{1}(2),srex_reg(1),'ver','top'   ,'hor','left' ,'color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_text(pol_lon{2}(2)   ,pol_lat{2}(2),srex_reg(2),'ver','top'   ,'hor','right','color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_text(pol_lon{3}(1)   ,pol_lat{3}(1),srex_reg(3),'ver','top'   ,'hor','left' ,'color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_text(pol_lon{4}(2)+15,pol_lat{4}(2),srex_reg(4),'ver','bottom','hor','left' ,'color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_text(pol_lon{5}(1)   ,pol_lat{5}(1),srex_reg(5),'ver','top'   ,'hor','left' ,'color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_text(pol_lon{6}(3)-10,pol_lat{6}(3),srex_reg(6),'ver','bottom','hor','right','color','r', 'Fontweight', 'Bold', 'Fontsize', 10); hold on;
m_line(lon_ispinb(ispinb_boundary), lat_ispinb(ispinb_boundary),'linewi',1,'color', colors(5,:)); % PINB area outline
export_fig figures/pct_irr_1915_lowres_noregions -transparent -png % save figure
% export_fig figures/pct_irr_1915 -m10 -transparent -png; close all; % save figure


% percentage of pixel equipped for irrigation map - year 2000
cbh = mf_plot_dom2(lon_mod, lat_mod, pct_irr, [], isunirr, caxes.Pirr, colormaps.Pirr, 0, 2, 'b', 'Area equipped for irrigation - 2000', 'Grid cell fraction equipped for irrigation [%]'); hold on;
set(cbh,'YTick',[0:10:50],'YTicklabel',[1 10:10:50])
export_fig figures/pct_irr_2000_lowres -transparent -png % save figure
% export_fig figures/pct_irr_2000 -m10 -transparent -png; close all; % save figure


% Population density
cbh = mf_plot_dom2(lon_mod, lat_mod, popdens, [], isinhabited, caxes.popd, colormaps.popd, 0, 2, 'c', 'Population density - 2000', 'Population density [inhabitants km^-^2]'); hold on;
set(cbh,'YTick',[0:100:500],'YTicklabel',[1 100:100:500])
export_fig figures/popdens_2000_lowres -transparent -png % save figure
% export_fig figures/popdens_2000 -m10 -transparent -png; close all; % save figure


% % Rural population density
% cbh = mf_plot_dom2(lon_mod, lat_mod, popdens_rural, [], isinhabited, caxes.popd, colormaps.popd, 0, 2, 'c', 'Rural population density - 2000', 'Rural population density [inhabitants km^-^2]'); hold on;
% set(cbh,'YTick',[0:100:500],'YTicklabel',[1 100:100:500])
% export_fig figures/popdens_rural_2000_lowres -transparent -png % save figure
% % export_fig figures/popdens_rural_2000 -m10 -transparent -png; close all; % save figure


end



% --------------------------------------------------------------------
% paper supplementary figure 2 - fig. 1 with older CRU version
% --------------------------------------------------------------------


if flags.plot_sfig2 == 1

    
disp('see python scripts from Auke Visser, reworked by Wim Thiery and Matthias Hauser');
disp('/net/exo/landclim/mathause/analysis/2018_for_wim_irrig/scripts/CRU_CESM_binnedboxplot_regional_WT_MH.ipynb');
disp('Mathias Huaser: I did it in a jupyter notebook - by executing ./startup it should launch the correct environment.');


end



% --------------------------------------------------------------------
% paper supplementary figure 3 - scatter plots with Delta TXm_irr removed
% --------------------------------------------------------------------


if flags.plot_sfig3 == 1

    
disp('see python scripts from Auke Visser, reworked by Wim Thiery and Matthias Hauser');
disp('/net/exo/landclim/mathause/analysis/2018_for_wim_irrig/scripts/CRU_CESM_binnedboxplot_regional_WT_MH.ipynb');
disp('Mathias Huaser: I did it in a jupyter notebook - by executing ./startup it should launch the correct environment.');


end



% --------------------------------------------------------------------
% paper supplementary figure 4 - irrigation-induced cooling from regression
% method
% --------------------------------------------------------------------


if flags.plot_sfig4 == 1



% CRU v4.02
mf_plot_dom2(lon_mod, lat_mod, CRUv402_t(:,:,1), [], island, caxes.regr, colormaps.regr, 0, 2, 'a', 'Observations (CRU v4.02)', 'Irrigation-induced change in TXm [K]'); hold on;
% export_fig figures/regr_CRUv402_lowres -transparent -png % save figure
export_fig figures/regr_CRUv402 -m10 -transparent -png; close all; % save figure


% CRU v3.22
mf_plot_dom2(lon_mod, lat_mod, CRUv322_t(:,:,1), [], island, caxes.regr, colormaps.regr, 0, 2, 'b', 'Observations (CRU v3.22)', 'Irrigation-induced change in TXm [K]'); hold on;
% export_fig figures/regr_CRUv322_lowres -transparent -png % save figure
export_fig figures/regr_CRUv322 -m10 -transparent -png; close all; % save figure
  
    
% CESM
mf_plot_dom2(lon_mod, lat_mod, CESM_t(:,:,1), [], island, caxes.regr, colormaps.regr, 0, 2, 'c', 'Model (CESM)', 'Irrigation-induced change in TXm [K]'); hold on;
% export_fig figures/regr_CESM_lowres -transparent -png % save figure
export_fig figures/regr_CESM -m10 -transparent -png; close all; % save figure
    

end



% --------------------------------------------------------------------
% paper supplementary figure 5 - risk ratio per srex region plots
% --------------------------------------------------------------------


if flags.plot_sfig5 == 1 && flags.refr == 0


% initialisation - y-axis limits (see srex_vars)
ylims = [1/6.25    6.25  ; ...  % WNA
         1/6.25    6.25  ; ...  % CNA
         1/6.25    6.25  ; ...  % MED
         1/6.25    6.25  ; ...  % WAS
         1/6.25    6.25  ; ...  % SAS
         1/6.25    6.25    ] ;  % EAS         

       
      
% prepare for loop
k = 1;


% loop oversrex regions
for i=1:6                      % all regions

    % construct matrices for plotting and rescale
    LR_bars=[LR_pic2ctl_srex_lp(:,i) LR_ctl2irr_srex_lp(:,i) LR_pic2irr_srex_lp(:,i)]; % all pixels!

    % determine whether legend needs to be plotted 
    if i == 1
        flag_plot_legend = 1;
    else
        flag_plot_legend = 0;
    end
    
    % plot srex region
    mf_plot_LR(LR_bars, percentages, ylims(i,:), alphabet(k), srex_reg{i}, flag_plot_legend)

    % save figure
    export_fig(sprintf('figures/LR_lp_%s_median.png',[num2str(k) '_' srex_reg{i}]), '-transparent');

    % for panel letters
    k = k + 1;
            
end




end



% --------------------------------------------------------------------
% paper supplementary figure 6 - binned changes
% --------------------------------------------------------------------


if flags.plot_sfig6 == 1

   
% initialisation - y-axis limits (see srex_vars)
ylims = [  0       0.33  ; ...  % WNA
           0       0.57  ; ...  % CNA
           0       0.92  ; ...  % MED
           0       0.64  ; ...  % WAS
           0       2.25  ; ...  % SAS
           0       0.21] ;      % EAS         
      

% prepare for loop
k = 1;


% bin other variables to TSA
for i=1                            %  QIRRIG only

    for j=1:6                      % MED only
        
        % get differences per bin    
        bindiff              = binmedian_binT_irr{i,j} - binmedian_binT_ctl{i,j};
        bindiff_Q25          = binQ25_binT_irr{i,j}    - binQ25_binT_ctl{i,j};
        bindiff_Q75          = binQ75_binT_irr{i,j}    - binQ75_binT_ctl{i,j};
        ind_pos              = find(bindiff > 0);
        ind_neg              = find(bindiff < 0);
        bindiff_pos          = zeros(size(bindiff));
        bindiff_neg          = zeros(size(bindiff));
        bindiff_pos(ind_pos) = bindiff(ind_pos);
        bindiff_neg(ind_neg) = bindiff(ind_neg);

        % quantile change of a given variable between CTL and IRR
        figure('OuterPosition',[100 200 950 410]);
        set(gcf, 'color', 'w');
        h(1) = bar(1:nbins,bindiff_pos,'histc'); hold on;
        h(2) = bar(1:nbins,bindiff_neg,'histc'); hold on;
        set(h(1),'facecolor',colors(17,:),'edgecolor','none');
        set(h(2),'facecolor',colors(17,:),'edgecolor','none');
        e = errorbar(1.5:nbins+0.5, bindiff, bindiff_Q25, bindiff_Q75,'k.','linewidth',1.5); hold on
        set(get(e,'children'),'clipping','off')
        axis([1 nbins+1 ylims(j,1) ylims(j,2)]);
        set(gca, 'Fontsize', 15, 'Fontweight', 'Bold');
        set(gca,'XTickLabel',{' ','P10',' ','P30',' ','P50',' ','P70', ' ','P90'}, 'Xtick', 1:nbins/10:nbins,'fontsize',15)
        xlabel('T_2_m percentiles', 'Fontsize', 18, 'Fontweight', 'Bold'); 
        ylabel(['\Delta ' srex_ylabels{i} ' [' srex_units{i} ']'], 'Fontsize', 18, 'Fontweight', 'Bold'); 
       
        text(1.15,ylims(j,2),alphabet(k),'ver','bottom','hor','center','Fontsize', 18)
        %text(nbins+1,ylims(i,2),[srex_vars{i} ', ' srex_reg{j}],'ver','bottom','hor','right','Fontsize', 11, 'Fontweight', 'Bold')
        text(nbins+1,ylims(j,2),srex_reg{j},'ver','bottom','hor','right','Fontsize', 18, 'Fontweight', 'Bold')

        % save figure
        export_fig(sprintf('figures/perc_change_%s.png',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');
        %export_fig(sprintf('figures/perc_change_%s.pdf',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');
                
        % for panel letters
        k = k + 1;
        
    end
    
end


end



% --------------------------------------------------------------------
% paper supplementary figure 7 - risk ratio box plots 
% --------------------------------------------------------------------


if flags.plot_sfig7 == 1 && flags.refr == 0

    
% initialisation - y-axis limits (see srex_vars)
ylims          = [1/6.25 6.25 ];
% ylims          = [1/64 16 ];


% prepare LR bar plots
LR_bars_lp   = [LR_pic2ctl_lp   LR_ctl2irr_lp   LR_pic2irr_lp  ];
LR_bars_ip   = [LR_pic2ctl_ip   LR_ctl2irr_ip   LR_pic2irr_ip  ];
LR_bars_pinb = [LR_pic2ctl_pinb LR_ctl2irr_pinb LR_pic2irr_pinb];


% plot all land
mf_plot_LR_boxplots(LR_bars_lp, boxplotdata_LR_pic2ctl_lp, boxplotdata_LR_ctl2irr_lp, boxplotdata_LR_pic2irr_lp, percentages, ylims, 'a', 'All land', 1)
export_fig('figures/LR_lp_boxplot.png', '-transparent');


% plot irrigated land
mf_plot_LR_boxplots(LR_bars_ip, boxplotdata_LR_pic2ctl_ip, boxplotdata_LR_ctl2irr_ip, boxplotdata_LR_pic2irr_ip, percentages, ylims, 'b', 'Irrigated land', 0)
export_fig('figures/LR_ip_boxplot.png', '-transparent');
  

% plot irrigated land
mf_plot_LR_boxplots(LR_bars_pinb, boxplotdata_LR_pic2ctl_pinb, boxplotdata_LR_ctl2irr_pinb, boxplotdata_LR_pic2irr_pinb, percentages, ylims, 'c', 'South Asia', 0)
export_fig('figures/LR_lp_PINB_boxplot.png', '-transparent');
      
    
    

end



% --------------------------------------------------------------------
% paper supplementary figure 8 - effect of irrigation in reference run
% --------------------------------------------------------------------


if flags.plot_sfig8 == 1

    
% prepare LR bar plots
LR_bars_pinb     = [LR_pic2ctl_pinb LR_ctl2irr_pinb LR_pic2irr_pinb];


% plot PINB - different depending on which reference runs
if     flags.refr == 0 % reference has irrigation (default)
    
    mf_plot_LR(LR_bars_pinb, percentages, [1/6.25 6.25], 'a', 'South Asia', 1)
    export_fig('figures/LR_median_pinb_sfig6_panela.png', '-transparent'); % save figure
    
elseif flags.refr == 1 % reference has no irrigation
    
    mf_plot_LR(LR_bars_pinb, percentages, [1/6.25 6.25], 'b', 'South Asia', 0)
    export_fig('figures/LR_median_pinb_sfig6_panelb.png', '-transparent'); % save figure
    
end


end



% --------------------------------------------------------------------
% paper supplementary figure 9 - change in hottest month
% --------------------------------------------------------------------


if flags.plot_sfig9 == 1

% get index of hottest month
[~, TREFHTMX_irr_hottestmonth] = nanmax(TREFHTMX_irr_mm, [], 3);
[~, TREFHTMX_pic_hottestmonth] = nanmax(TREFHTMX_pic_mm, [], 3);

% Hottest month - 20cirr
cbh = mf_plot_dom2(lon_mod, lat_mod, TREFHTMX_pic_hottestmonth, [], island, caxes.mon, colormaps.mon, 0, 2, 'a', 'Hottest month (1901-1930)', 'Hottest month'); hold on;
cbh.Ticks      = 1:12;
cbh.TickLabels = {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'};
export_fig figures/hottest_month_20cirr_lowres -transparent -png % save figure
% export_fig figures/hottest_month_20cirr -m10 -transparent -png; close all; % save figure

% Hottest month - change
mf_plot_dom2(lon_mod, lat_mod, TREFHTMX_irr_hottestmonth - TREFHTMX_pic_hottestmonth, [], island, caxes.dmon, colormaps.dmon, 0, 2, 'b', 'Change in hottest month (1901-1930 to 1981-2010)', 'Change in hottest month'); hold on;
export_fig figures/hottest_month_change_lowres -transparent -png % save figure
% export_fig figures/hottest_month_change -m10 -transparent -png; close all; % save figure

    
end



% --------------------------------------------------------------------
% survival plots
% Note: commented out "fix_lines(name);" in print2eps    
% --------------------------------------------------------------------


if flags.plot_surv == 1
        

% set survival plots axes limits
surv_ylims = [0.00002 0.1]; %[0.00002 0.5]


% initialisation - y-axis limits (see srex_vars)
surv_xlims = [299     305     ; ...  % QIRRIG
              299     305     ; ...  % TSA
              299     305     ; ...  % Qle
              299     305     ; ...  % FSH
              299     305     ; ...  % FSA
              299     305     ; ...  % FIRA
              299     305     ; ...  % PRECT
              299     305     ; ...  % PRECT
              299     305     ; ...  % PRECT
              299     305     ; ...  % PRECT
              299     305   ] ;      % QRUNOFF         
      

% prepare for loop
k = 1;


% loop over srex variables
% for i=[2 9 10 11]                 % all variables
for i=11                        % TREFHTMX only

    % loop over srex regions
%     for j=1:7                  % all regions    
    for j=5                      % SAS only


        % survival plots of daily accumulated precipitation - Lake Victoria
        figure%('OuterPosition',[150 200 300 450]);
        set(gcf, 'color', 'w');
        set(gca,'color','w');
        semilogy(surv_x_pic{i,j}, surv_f_pic{i,j}, 'k'  , 'LineWidth',1.5); hold on;
        semilogy(surv_x_ctl{i,j}, surv_f_ctl{i,j}, 'r'  , 'LineWidth',1.5); hold on;
        semilogy(surv_x_irr{i,j}, surv_f_irr{i,j}, 'b'  , 'LineWidth',1.5); hold on;
        set(gca, 'Fontsize', 12, 'Fontweight', 'Bold');
        set(gca, 'ylim', surv_ylims);
%         axis([surv_xlims(i,1) surv_xlims(i,2) surv_ylims(1) surv_ylims(2)]);        
        xlabel([srex_ylabels{i} ' [' srex_units{i} ']'], 'Fontsize', 18, 'Fontweight', 'Bold'); 
        ylabel('Normalized frequency', 'Fontweight', 'Bold', 'Fontsize', 12);
        % set(gca, 'XTick', [1, 10, 20, 30, 40, 50 60])
        % set(gca, 'YTick', [0.02 0.03 0.05 0.1 0.2 0.3 0.5 1.0])
%         legend('CTL','IRR', 3);
        legend('PIC','CTL','IRR', 3);
        set(legend, 'YColor', 'w', 'XColor', 'w', 'Fontweight', 'Bold', 'Fontsize', 10);
        text(1,1.22,{'a'},'ver','top','hor','left', 'Fontweight', 'Bold', 'Fontsize', 11); 
        text(60,1.22,{'SAS'},'ver','top','hor','right', 'Fontweight', 'Bold', 'Fontsize', 11); 
        hold off;
        box on;


        % plot panel letter and region
        xlims_plot = get(gca, 'xlim');
        ylims_plot = get(gca, 'ylim');
        text(xlims_plot(1),ylims_plot(2),alphabet(k),'ver','bottom','hor','left' ,'Fontsize', 18)
        text(xlims_plot(2),ylims_plot(2),srex_reg{j},'ver','bottom','hor','right','Fontsize', 18, 'Fontweight', 'Bold')


        
        % save figure
        export_fig(sprintf('figures/surv_c_%s.png',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');
        %export_fig(sprintf('figures/surv_%s.pdf',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');

               
        % for panel letters
        k = k + 1;
        
    end
    
end


end



% --------------------------------------------------------------------
% risk ratio plots  
% --------------------------------------------------------------------


if flags.plot_LRts == 1
            


% initialisation - y-axis limits (see srex_vars)
ylims = [ 1/8      8     ; ...  % WNA
          1/8      8     ; ...  % CNA
          1/8      8     ; ...  % MED
          1/8      8     ; ...  % WAS
          1/8      8     ; ...  % SAS
          1/8      8     ; ...  % SEA
          1/8      8       ] ;  % EAS         

       
% get x-axis ticklabels
for i=1:npercentages
    xticklabels{i} = ['>' num2str(percentages(i)) '%']; %#ok<SAGROW>
end


% prepare for loop
k = 1;


% bin other variables to TSA
% for i=1:11                        % all variables
for i=11                        % TREFHTMX only

    for j=1:7                      % all regions
%     for j=5                      % SAS only
        
        
        % construct matrices for plotting and rescale
        [LR_ts_plot       , LR_ts_plot_pos, LR_ts_plot_neg] = mf_riskratio_plotscale([LR_ts_pic2ctl{i,j}; LR_ts_ctl2irr{i,j}; LR_ts_pic2irr{i,j}]);
%         [LR_ts_plot_ci_low, ~             , ~             ] = mf_riskratio_plotscale([LR_ts_ci_pic2ctl{i,j}(:,1) LR_ts_ci_ctl2irr{i,j}(:,1) LR_ts_ci_pic2irr{i,j}(:,1)]');
%         [LR_ts_plot_ci_up , ~             , ~             ] = mf_riskratio_plotscale([LR_ts_ci_pic2ctl{i,j}(:,2) LR_ts_ci_ctl2irr{i,j}(:,2) LR_ts_ci_pic2irr{i,j}(:,2)]');
        

        % quantile change of a given variable between CTL and IRR
        figure%('OuterPosition',[100 200 950 410]);
        set(gcf, 'color', 'w');
        h = bar(1:npercentages, LR_ts_plot'); hold on;
        set(h(1),'BaseValue',1)
        set(h,'edgecolor','none');
        colormap([colors(16,:); colors(17,:); colors(15,:)])
%         d = errorbar(1+0.1:npercentages+0.1, LR_ts_plot(1,:), LR_ts_plot_ci_low(1,:), LR_ts_plot_ci_low(1,:),'k.','linewidth',1.5); hold on
%         e = errorbar(1+0.3:npercentages+0.3, LR_ts_plot(2,:), LR_ts_plot_ci_low(2,:), LR_ts_plot_ci_low(2,:),'k.','linewidth',1.5); hold on
%         f = errorbar(1+0.55:npercentages+0.55, LR_ts_plot(3,:), LR_ts_plot_ci_low(3,:), LR_ts_plot_ci_low(3,:),'k.','linewidth',1.5); hold on
%         set(get(d,'children'),'clipping','off')
%         set(get(e,'children'),'clipping','off')
%         set(get(f,'children'),'clipping','off')
        axis([0.5 npercentages+0.5 mf_riskratio_plotscale(ylims(j,1)) mf_riskratio_plotscale(ylims(j,2))]);
        set(gca, 'Fontsize', 13, 'Fontweight', 'Bold', 'Xcolor', axcolor, 'Ycolor', axcolor);
        mf_xticklabel_rotate(1.0:npercentages+0.0,90,xticklabels,'interpreter','none','Fontsize', 11, 'Fontweight', 'Bold', 'color', axcolor); 
        set(gca,'YTickLabel',{'1/8','1/6','1/4','1/2','1','2','4','6','8'}, 'Ytick', mf_riskratio_plotscale([1/8 1/6 1/4 1/2 1 2 4 6 8]))
        xlabel('T_{2m,max} percentiles', 'Fontsize', 18, 'Fontweight', 'Bold'); 
        ylabel('LR', 'Fontsize', 18, 'Fontweight', 'Bold'); 
        legend('present-day warming','present-day irrigation','present-day warming + irrigation', 2);
        set(legend, 'YColor', 'w', 'XColor', 'w', 'Fontweight', 'Bold', 'Fontsize', 10, 'textcolor', axcolor);
        text(0.6,ylims(j,2),alphabet(k),'ver','bottom','hor','center','Fontsize', 18, 'color', axcolor)
        text(npercentages+0.5,ylims(j,2),srex_reg{j},'ver','bottom','hor','right','Fontsize', 18, 'Fontweight', 'Bold', 'color', axcolor)


        % save figure
        export_fig(sprintf('figures/LR_ts_%s_pic150yrs.png',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');
        %export_fig(sprintf('figures/LR_ts_%s.pdf',[num2str(k) '_' srex_vars{i} '_' srex_reg{j}]), '-transparent');
                
        % for panel letters
        k = k + 1;
        
    end
    
end

end




