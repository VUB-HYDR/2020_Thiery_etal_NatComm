% --------------------------------------------------------------------
% function to select all pixels within a certain SREX region and
% then compute the spatial mean using the mf_fieldmean function
% SREX regions are defined in Seneviratne et al. (2012: fig 3.1)
% --------------------------------------------------------------------


function [out, out_pol_lat, out_pol_lon] = mf_srex(lat, lon, var, user_mask, varargin)

% input to the function is an input variable (2D or 3D) and its lat/lon  
% coordinates, and then 'region label' as defined by the SREX report.
% it is possible to input several srex regions in the same function call

% this function assumes that the input grid is 2D, and spanning -90N to +90N 
% and -180E to 180 E. It also assumes a 'global' land mask called 'island'

% % TO DO
% the function outputs NaN if data is missing more than 50% of the land
% pixels (CESM grid as reference)

% e.g. to get the average 2m temperature for East and South Africa
% [T2M_EAF, T2M_SAF] = mf_srex(lat_mod, lon_mod, T_2M, 'EAF', 'SAF');


% invoke globals
global island



% --------------------------------------------------------------------
% Initialisation
% --------------------------------------------------------------------


% Define srex region labels
labels = {'ALA', 'AMZ', 'CAM', 'CAS', 'CEU', 'CGI', 'CNA', 'EAF', 'EAS', 'ENA', 'MED', 'NAS', 'NAU', 'NEB', 'NEU', 'SAF', 'SAH', 'SAS', 'SAU', 'SSA', 'SEA', 'TIB', 'WAF', 'WAS', 'WSA', 'WNA'};


% Define coordinates (Latitude [°], Longitude [°]) of Region Corners
% see Table 3.A-1 (Appendix 3.A: Notes and Technical Details on Chapter 3 Figures)
coordinates = {...
'(60.000N, 105.000W)', '(60.000N, 168.022W)', '(72.554N, 168.022W)', '(72.554N, 105.000W)' , ''                   , ''                  ; ...
'(20.000S, 66.377W)' , '(1.239S, 79.729W)'  , '(11.439N, 68.800W)' , '(11.439N, 50.000W)'  , '(20.000S, 50.000W)' , ''                  ; ...
'(11.439N, 68.800W)' , '(1.239S, 79.729W)'  , '(28.566N, 118.323W)', '(28.566N, 90.315W)'  , ''                   , ''                  ; ...
'(30.000N, 60.000E)' , '(50.000N, 60.000E)' , '(50.000N, 75.000E)' , '(30.000N, 75.000E)'  , ''                   , ''                  ; ...
'(45.000N, 10.000W)' , '(48.000N, 10.000W)' , '(61.320N, 40.000E)' , '(45.000N, 40.000E)'  , ''                   , ''                  ; ...
'(50.000N, 10.000W)' , '(50.000N, 105.000W)', '(85.000N, 105.000W)', '(85.000N, 10.000W)'  , ''                   , ''                  ; ...
'(50.000N, 85.000W)' , '(28.566N, 85.000W)' , '(28.566N, 105.000W)', '(50.000N, 105.000W)' , ''                   , ''                  ; ...
'(11.365S, 25.000E)' , '(15.000N, 25.000E)' , '(15.000N, 51.990E)' , '(11.365S, 51.990E)'  , ''                   , ''                  ; ...
'(20.000N, 100.000E)', '(50.000N, 100.000E)', '(50.000N, 145.000E)', '(20.000N, 145.000E)' , ''                   , ''                  ; ...
'(25.000N, 60.000W)' , '(25.000N, 85.000W)' , '(50.000N, 85.000W)' , '(50.000N, 60.000W)'  , ''                   , ''                  ; ...
'(30.000N, 10.000W)' , '(45.000N, 10.000W)' , '(45.000N, 40.000E)' , '(30.000N, 40.000E)'  , ''                   , ''                  ; ...
'(50.000N, 40.000E)' , '(70.000N, 40.000E)' , '(70.000N, 180.000E)', '(50.000N, 180.000E)' , ''                   , ''                  ; ...
'(30.000S, 110.000E)', '(10.000S, 110.000E)', '(10.000S, 155.000E)', '(30.000S, 155.000E)' , ''                   , ''                  ; ...
'(20.000S, 34.000W)' , '(20.000S, 50.000W)' , '(0.000N, 50.000W)'  , '(0.000N, 34.000W)'   , ''                   , ''                  ; ...
'(48.000N, 10.000W)' , '(75.000N, 10.000W)' , '(75.000N, 40.000E)' , '(61.320N, 40.000E)'  , ''                   , ''                  ; ...
'(35.000S, 10.000W)' , '(11.365S, 10.000W)' , '(11.365S, 51.990E)' , '(35.000S, 51.990E)'  , ''                   , ''                  ; ...
'(15.000N, 20.000W)' , '(30.000N, 20.000W)' , '(30.000N, 40.000E)' , '(15.000N, 40.000E)'  , ''                   , ''                  ; ...
'(5.000N, 60.000E)'  , '(30.000N, 60.000E)' , '(30.000N, 100.000E)', '(20.000N, 100.000E)' , '(20.000N, 95.000E)' , '(5.000N, 95.000E)' ; ...
'(50.000S, 110.000E)', '(30.000S, 110.000E)', '(30.000S, 180.000E)', '(50.000S, 180.000E)' , ''                   , ''                  ; ...
'(20.000S, 39.376W)' , '(56.704S, 39.376W)' , '(56.704S, 67.348W)' , '(50.000S, 72.141W)'  , '(20.000S, 66.377W)' , ''                  ; ...
'(10.000S, 95.000E)' , '(20.000N, 95.000E)' , '(20.000N, 155.000E)', '(10.000S, 155.000E)' , ''                   , ''                  ; ...
'(30.000N, 75.000E)' , '(50.000N, 75.000E)' , '(50.000N, 100.000E)', '(30.000N, 100.000E)' , ''                   , ''                  ; ...
'(11.365S, 20.000W)' , '(15.000N, 20.000W)' , '(15.000N, 25.000E)' , '(11.365S, 25.000E)'  , ''                   , ''                  ; ...
'(15.000N, 40.000E)' , '(50.000N, 40.000E)' , '(50.000N, 60.000E)' , '(15.000N, 60.000E)'  , ''                   , ''                  ; ...
'(1.239S, 79.729W)'  , '(20.000S, 66.377W)' , '(50.000S, 72.141W)' , '(56.704S, 67.348W)'  , '(56.704S, 82.022W)' , '(0.530N, 82.022W)' ; ...
'(28.566N, 105.000W)', '(28.566N, 130.000W)', '(60.000N, 130.000W)', '(60.000N, 105.000W)' , ''                   , ''                       };


% Minimum number of land pixels that must contain a value in order to be
% retained
min_nobs = 50; % [%]


% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% check if user_mask is defined
if isempty(user_mask)
    user_mask = island; % safety
end


% prepare for loop
out = NaN(length(varargin),1);


% loop over input srex regions
for i=1:length(varargin)

    
    % set the region
    label = varargin{i};
    
    
    % get index corresponding to the region
    ind = strfind(labels, label);
    ind = find(not(cellfun('isempty', ind)));

    
    % get coordinates of corners
    corners  = coordinates(ind,:);                                          %#ok<FNDSB>
    ncorners = find(not(cellfun('isempty', corners)), 1, 'last'); % number of corners 
    corners  = corners(1:ncorners);                               % remove empty cells
    
    
    % get polygon of the region
    pol_lat = NaN(ncorners, 1);
    pol_lon = NaN(ncorners, 1);
    for j=1:ncorners
                
        % get lat/lon coordinates of each corner
        ind_sep = strfind(corners{1},','); % index of the comma separator
        pol_lat(j) = str2double(corners{j}(2         : ind_sep-2));
        pol_lon(j) = str2double(corners{j}(ind_sep+2 : end    -2));
        
        % invert lat coordinate sign if corner is in the southern hemisphere
        if corners{j}(ind_sep-1) == 'S'
            pol_lat(j) = pol_lat(j) .* -1;
        end
            
        % invert lon coordinate sign if corner is in the western hemisphere
        if corners{j}(end-1) == 'W'
            pol_lon(j) = pol_lon(j) .* -1;
        end
        
    end
    
    
    % select pixels inside polygon
    mask = inpolygon(lon, lat, pol_lon, pol_lat);

    
    % retain land pixels only and apply user-defined mask
    mask = mask & island & user_mask;
    
    
    % % debugging: plot mask
    % mf_plot_dom2(lon, lat, double(mask), [0 1], mf_colormap_cpt('BrBG_03', 2), 1.9, 0, 1, ' ', label, ' ');
 

    % get mean over the srex region (land pixels only)
    [~, var_srex] = mf_fieldmean(var, mask);

    
    % in case var has a third dimension (e.g. time), get mean over that dimension as well
    var_srex = nanmean(var_srex);
    
    
    % retain value only if data was available for more than X% of all land pixels in the subdomain 
    npixels   = length(find(       var(mask)));
    nobs      = length(find(~isnan(var(mask))));
    nobs_perc = (nobs / npixels) * 100;
    if nobs_perc < min_nobs;
        var_srex = NaN;
    end
    
    
    % store data
    % varargout{i} = var_srex; % if you use varargout in function call
    out(i) = var_srex;
    out_pol_lat{i} = pol_lat; %#ok<*AGROW>
    out_pol_lon{i} = pol_lon;

    
end


end

