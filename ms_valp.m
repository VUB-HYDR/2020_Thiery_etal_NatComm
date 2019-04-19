
% --------------------------------------------------------------------
% Compute values used in the paper
% note: preferably run "main"
% --------------------------------------------------------------------


clc;


% --------------------------------------------------------------------
% Results
% --------------------------------------------------------------------


% percentage of land area with increase in irrigated fraction above 35%
isirr_gt35perc       = (pct_irr - pct_irr_1915) > 35;
area_irr_gt35perc    = nansum(nansum(area .* isirr_gt35perc)); % [km^2]
area_lp              = nansum(nansum(area .* island        )); % [km^2]
percent_irr_gt35perc = area_irr_gt35perc ./ area_lp .* 100


% regional median LR from global warming
for i=1:nreg
    
    % get ndays expected per year in early 20h century (trivial)
    [~, range_ndays_srex_P99_ref_lp(:,i)    , range_ndays_srex_P99_ref_ip(:,i)    ] = mf_fieldmedian(prob_ref(:,:,5)     .* 365, issrex(:,:,i) & island, issrex(:,:,i) & isirr); %#ok<SAGROW>
    
    % get ndays expected per year today
    [~, range_ndays_srex_P99_pic2ctl_lp(:,i), range_ndays_srex_P99_pic2ctl_ip(:,i)] = mf_fieldmedian(prob_pic2ctl(:,:,5) .* 365, issrex(:,:,i) & island, issrex(:,:,i) & isirr); %#ok<SAGROW>

end
range_LR_srex_P99 = LR_pic2ctl_srex_lp(5,:)
range_ndays_srex_P99_ref_lp
range_ndays_srex_P99_pic2ctl_lp



% --------------------------------------------------------------------
% Discussion
% --------------------------------------------------------------------


% number of people half as much exposed to higher hot extremes (TX99p) in the year 2000 thanks to irrigation
% first compute the exosure mask: 
isexposed                = LR_ctl2irr(:,:,5) < 0.5;
isexposed(lat_mod < -60) = 0;                       % remove antarctica
isexposed(~island)       = 0;                       % remove oceans


% then compute global population and number of people exposed
npeople_ap               = nansum(nansum(             popdens .* area .* pct_land./100)); % [#people/km^2 * km^2 = #people]
npeople_exposed          = nansum(nansum(isexposed .* popdens .* area .* pct_land./100))  % [#people/km^2 * km^2 = #people]


% then compute global RURAL population and number of RURAL people exposed
npeople_rural_ap         = nansum(nansum(             popdens_rural .* area .* pct_land./100)); % [#people/km^2 * km^2 = #people]
npeople_rural_exposed    = nansum(nansum(isexposed .* popdens_rural .* area .* pct_land./100))  % [#people/km^2 * km^2 = #people]



% --------------------------------------------------------------------
% Methods
% --------------------------------------------------------------------


% compute irrigated area in 1915 and 200
area_irr_1915 = nansum(nansum(area .* pct_irr_1915./100)) % [km^2]
area_irr_2000 = nansum(nansum(area .* pct_irr     ./100)) % [km^2]




% --------------------------------------------------------------------
% numbers for review
% --------------------------------------------------------------------


% number of pixels with irrigated crop fraction above 70%  !! at 1 degree but fig 1 is at 0.5 degree (but using second-order conservative remapping) !!
npixels_pct_irr_gt70 = length(pct_irr(pct_irr > 70));


