
% --------------------------------------------------------------------
% subroutine to load the variables
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% Load observational data
% --------------------------------------------------------------------


% load UNWPP population density data (GPW v4) for the year 2000 [inh/km^2]
[~, ~, popdens] = mf_load('UNWPP_v4_2000-2010_yearmean.nc' , 'UN-Adjusted Population Density, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes');
popdens         = popdens(:,:,1); % select year 2000


% load rural popultion density from HYDE3.2 for the year 2000 [inh/km^2]
[~, ~, popdens_rural] = mf_load('HYDE32_rpopd_2000_yearmean.nc', 'rurc_');



% --------------------------------------------------------------------
% Load CESM model data
% --------------------------------------------------------------------


% load 2D model constants (pct_irr = percent of grid cell that is irrigated, see mksrfdat.F90)
[lat_mod, lon_mod, pct_irr, pct_land, area, pct_pft] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.001_constants.nc', 'PCT_IRRIG', 'landfrac', 'AREA', 'PCT_PFT');
pct_irr_1915                                         = circshift(rot90(ncread('surfdata_0.9x1.25_simyr2000_irrcr_c100916_irr1915.nc', 'PCT_IRRIG')), [0 size(lat_mod,2)./2]);


% load 2D model atm variables
[~, ~, TREFHT_ctl, TXx_ctl, TNn_ctl, CDD_ctl, HWDI_ctl, WSDI_ctl] = mf_load('f.e122.F1850PDC5.f09_g16.control-io192.ensmean_atm_h1_timmean.nc'   , 'TREFHT', 'TXx', 'TNn', 'CDD', 'HWDI', 'WSDI');
[~, ~, TREFHT_irr, TXx_irr, TNn_irr, CDD_irr, HWDI_irr, WSDI_irr] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.ensmean_atm_h1_timmean.nc', 'TREFHT', 'TXx', 'TNn', 'CDD', 'HWDI', 'WSDI');


% load 2D model lnd variables
[~, ~, TSA_ctl, Rx1day_ctl, Rx5day_ctl] = mf_load('f.e122.F1850PDC5.f09_g16.control-io192.ensmean_lnd_h1_timmean.nc'   , 'TSA', 'Rx1day', 'Rx5day');
[~, ~, TSA_irr, Rx1day_irr, Rx5day_irr] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.ensmean_lnd_h1_timmean.nc', 'TSA', 'Rx1day', 'Rx5day');


% load 2D model lnd variables
[~, ~, TREFHTMX_irr_mm] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.ensmean_atm_h1_ymonmean.nc', 'TREFHTMX');
[~, ~, TREFHTMX_pic_mm] = mf_load('f.e122.F1850PDC5.f09_g16.20cirr-io192.ensmean_atm_h1_ymonmean.nc'    , 'TREFHTMX');


% load data averaged over srex regions for survival plots
if flags.surv == 1  
for i=1:length(srex_vars)    
    for j=1:length(srex_reg)
        
        % Pre-Industrial Control (CMIP-type)
        try
            srex_pic{i,j} = squeeze(ncread(['b.e122.B1850C5CN.f19_g16.control-io144.500_' srex_block{i} '_h1_srex_land.nc'], [srex_vars{i} '_' srex_reg{j}])); %#ok<*SAGROW>
        catch errmessage
            srex_pic{i,j} = 0;
        end
        
        % Present-day (AMIP-type)
        for k=1:nens
            srex_ctl{i,j,k} = squeeze(ncread(['f.e122.F1850PDC5.f09_g16.control-io192.00'    num2str(k) '_' srex_block{i} '_h1_srex_land.nc'], [srex_vars{i} '_' srex_reg{j}])); %#ok<*SAGROW>
            srex_irr{i,j,k} = squeeze(ncread(['f.e122.F1850PDC5.f09_g16.irrigation-io192.00' num2str(k) '_' srex_block{i} '_h1_srex_land.nc'], [srex_vars{i} '_' srex_reg{j}]));
        end

    end
end
end


% Load percentage maps 
for i=1:length(percentages)    

    % reference maps
    perc_ref(:,:,i) = ones(size(lat_mod)) .* percentages(i);
    

    % Pre-Industrial Control (pic) to present-day control (ctl)
    if     flags.refr == 0 % reference has irrigation (default)
        [~, ~, perc_pic2ctl(:,:,i)] = mf_load('f.e122.F1850PDC5.f09_g16.control-io192.merged_atm_h1_TREFHTMX_prob_20cirr2ctl.nc'   , ['prob_pic2ctl_P' num2str(percentages(i),'%.3f')]);
    elseif flags.refr == 1 % reference has no irrigation
        [~, ~, perc_pic2ctl(:,:,i)] = mf_load('f.e122.F1850PDC5.f09_g16.control-io192.merged_atm_h1_TREFHTMX_prob_20cc2ctl.nc'   , ['prob_pic2ctl_P' num2str(percentages(i),'%.3f')]);
    end
       
    
    % Present-day control (ctl) to present-day irrigation (irr)
    [~, ~, perc_ctl2irr(:,:,i)] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.merged_atm_h1_TREFHTMX_prob_ctl2irr.nc', ['prob_ctl2irr_P' num2str(percentages(i),'%.3f')]);

           
    % Pre-Industrial Control (pic) to present-day irrigation (irr)
    if     flags.refr == 0 % reference has irrigation (default)
        [~, ~, perc_pic2irr(:,:,i)] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.merged_atm_h1_TREFHTMX_prob_20cirr2irr.nc', ['prob_pic2irr_P' num2str(percentages(i),'%.3f')]);
    elseif flags.refr == 1 % reference has no irrigation
        [~, ~, perc_pic2irr(:,:,i)] = mf_load('f.e122.F1850PDC5.f09_g16.irrigation-io192.merged_atm_h1_TREFHTMX_prob_20cc2irr.nc', ['prob_pic2irr_P' num2str(percentages(i),'%.3f')]);
    end

end


% load spatially-averaged lnd variables for binning
if flags.perc == 1  
for i=1:length(srex_vars)    
    for j=1:length(srex_reg)
        for k=1:nens
            srex_ctl{i,j,k} = squeeze(ncread(['f.e122.F1850PDC5.f09_g16.control-io192.00'    num2str(k) '_' srex_block{i} '_h1_srex.nc'], [srex_vars{i} '_' srex_reg{j}])); %#ok<*SAGROW>
            srex_irr{i,j,k} = squeeze(ncread(['f.e122.F1850PDC5.f09_g16.irrigation-io192.00' num2str(k) '_' srex_block{i} '_h1_srex.nc'], [srex_vars{i} '_' srex_reg{j}]));
        end
    end
end
end



% --------------------------------------------------------------------
% Load results from observational analysis
% --------------------------------------------------------------------


% irrigation-induced cooling according to regression method (Lejeune et al., 2018)
% first dimension represent irrigation thresholds 0, 0.1, ..., 0.8

% load CRU_t v3.22
load('formatlab_CRUv322_t.mat');  
CRUv322_t  = flipud(permute(CRU_t, [2 3 1]));


% load CRU_t v4.02
load('formatlab_CRUv402_t.mat'); 
CRUv402_t  = flipud(permute(CRU_t, [2 3 1]));


% load CESM_t
load('formatlab_CESM_t.mat'); 
CESM_t     = flipud(permute(CESM_t   , [2 3 1]));



% --------------------------------------------------------------------
% unit conversions
% --------------------------------------------------------------------


% pct_land
pct_land = pct_land .* 100; % [] to [%]
