

% --------------------------------------------------------------------
% function to plot domain data
% --------------------------------------------------------------------


function [cbh] = mf_plot_dom2(lon, lat, var, var_sign, island, caxis_dom, colormap_dom, flag_sp, flag_cb, panel, experiment, cbcaption)


% note
% flag_cb = 0; no colorbar
%           1: standard colorbar
%           2: thin colorbar southoutside



% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% Set grid characteristics - whole domain
lat_min_dom =  -70.0;
lat_max_dom =   80.0;
lon_min_dom = -177.5;
lon_max_dom =  177.5;


% set model resolution
res_lat = 0.9; 
res_lon = 1.25;


% set sea color
seacolor = [0.95 0.95 0.95]; % very light grey


% define axes color                               
axcolor = [0.3 0.3 0.3]; % 70% contrast (so 0.3) is advised



% --------------------------------------------------------------------
% manipulations: hatching
% hatch insignificant pixels only if the sum of their indices is even
% --------------------------------------------------------------------


% if no significance mask is given, assume all pixels significant
if isempty(var_sign)
    var_sign = ones(size(var));
end


% retain only insignificant values on land
var_insign                = NaN(size(var_sign));
var_insign(var_sign == 0) = 1;


% do not consider potential significant changes over oceans
var_insign(~island)       = NaN; 


% generate checkerboard matrix
[nlat, nlon]                     = size(lat);
checkerboard                     = zeros(nlat, nlon);   % "wide" checkerboard - perfect
checkerboard(1:3:nlat, 1:3:nlon) = 1;
checkerboard(3:3:nlat, 2:3:nlon) = 1;
checkerboard(2:3:nlat, 3:3:nlon) = 1;


% to have less lines, hatch insignificant pixels only if the sum of their indices is even
var_insign(checkerboard == 0) = NaN; 


% prepare for hatching out
lon_insign = lon.*var_insign;
lat_insign = lat.*var_insign;
lon_insign = lon_insign(~isnan(lon_insign));
lat_insign = lat_insign(~isnan(lat_insign));
lon_insign = [lon_insign' - res_lon/2; lon_insign' + res_lon/2]; % plotting these gives diagonal over the pixel
lat_insign = [lat_insign' - res_lat/2; lat_insign' + res_lat/2];



% --------------------------------------------------------------------
% Visualisation
% --------------------------------------------------------------------


% design figure
if flag_sp == 0
figure('OuterPosition',[100 200 800 500]);
set(gcf, 'color', 'w');
set(gca,'color','w');
end


% Initialise grid and projection
m_proj('Equidistant Cylindrical','long',[lon_min_dom lon_max_dom],'lat',[lat_min_dom lat_max_dom]);                                        % generate Equidistant Cylindrical projection (alternatives: Miller or Mercator)
m_grid('box', 'on', 'xtick', [], 'ytick', [], 'color', axcolor, 'Fontweight', 'Bold', 'FontSize', 11, 'linewidth', 0.5, 'linestyle', 'none'); hold on; % draw grid with normal box


% set shading flat and ignore warnings this generates - necessary only for writing text on maps
shading flat;                                                                                                      % set shading flat so that you can still plot text on the map                                          
warning('off','all');


% plot light grey background
var(island==0) = NaN;
hback          = m_pcolor(lon - res_lon/2, lat + res_lat/2, ones(size(var)));  hold on;         
set(hback, 'edgecolor', 'none', 'facecolor', seacolor); 


% % plot country borders 
% M = m_shaperead('ne_10m_admin_0_countries');                                                                       % load country borders
% for k=1:length(M.ncst), 
%   m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [0.6 0.6 0.6]); hold on;                                          % plot country borders
% end;
% hold on;


% plot data
g = m_pcolor(lon - res_lon/2, lat + res_lat/2, var);                                                               % with shading flat will draw a panel between the (i,j),(i+1,j),(i+1,j+1),(i,j+1) coordinates of the lon/lat matrices with a color corresponding to the data value at (i,j).
set(g, 'edgecolor', 'none');                                                                                       % remove black grid around pixels 
set(gca, 'Fontsize', 14, 'Fontweight', 'Bold');                                                                    % set axes properties
caxis(caxis_dom);                                                                                                  % set colorscale axes
colormap(colormap_dom); mf_freezeColors;                                                                           % select colormap and freeze it
if flag_cb == 1
    cbh=colorbar('Fontsize', 14, 'Fontweight', 'Bold', 'color', axcolor); hold on;     
    cbfreeze; mf_freezeColors;
elseif flag_cb == 2
    cbh=colorbar('location', 'SouthOutside', 'Fontsize', 14, 'Fontweight', 'Bold', 'TickLength', [0 0], 'color', axcolor); hold on;     
    y1=get(gca,'position');
    y=get(cbh,'Position');
    y(4)=0.03;
    set(cbh,'Position',y)
    set(gca,'position',y1)
    mf_freezeColors;
end
hold on;


% plot hatching
m_plot(lon_insign,lat_insign, '-', 'color', [0.8 0.8 0.8], 'linewidth', 1); hold on;                            % hatch over insignificant areas


% plot coasts
m_gshhs('cc','color', axcolor);                                                                                         % add coastline


% add text
m_text(-175,81,panel,'ver','bottom','hor','left', 'Fontweight', 'Bold', 'Fontsize', 12, 'color', axcolor); hold on;
m_text(175,81,experiment,'ver','bottom','hor','right', 'Fontweight', 'Bold', 'Fontsize', 11, 'color', axcolor); hold on;
if flag_cb == 1
    m_text(240,0,cbcaption,'ver','bottom','hor','center', 'Fontweight', 'Bold', 'Fontsize', 12,'rotation',-90, 'color', axcolor); hold off;
elseif flag_cb == 2
    m_text(0,-135,cbcaption,'ver','bottom','hor','center', 'Fontweight', 'Bold', 'Fontsize', 12, 'color', axcolor); hold off;
end


end
