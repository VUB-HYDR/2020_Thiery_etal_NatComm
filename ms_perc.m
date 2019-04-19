
% --------------------------------------------------------------------
% Percentile change routine
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% initialise bin characteristics
nbins = 20;



% --------------------------------------------------------------------
% manipulations: binned changes
% --------------------------------------------------------------------


% binning
for i=1:length(srex_vars)    
    for j=1:nreg

        
        % bin all variables according to PRECT        
        for k=1:nens
            
            % get TSA bins (use "tabulate(bin)" for some bin information)
            binT_ctl = ceil(nbins * tiedrank(srex_ctl{2,j,k}) / length(srex_ctl{2,j,k})); % bin according to TSA (temperature)
            binT_irr = ceil(nbins * tiedrank(srex_irr{2,j,k}) / length(srex_irr{2,j,k})); % bin according to TSA (temperature)
            
            % binning
            [dummy_mean_ctl(:,k), dummy_median_ctl(:,k), ~, dummy_Q25_ctl(:,k), dummy_Q75_ctl(:,k)] = mf_bin(srex_ctl{i,j,k}, binT_ctl, nbins); %#ok<*SAGROW>
            [dummy_mean_irr(:,k), dummy_median_irr(:,k), ~, dummy_Q25_irr(:,k), dummy_Q75_irr(:,k)] = mf_bin(srex_irr{i,j,k}, binT_irr, nbins);
            
        end        
        binmean_binT_ctl{i,j}   = nanmean(dummy_mean_ctl  , 2);
        binmean_binT_irr{i,j}   = nanmean(dummy_mean_irr  , 2);
        binmedian_binT_ctl{i,j} = nanmean(dummy_median_ctl, 2);
        binmedian_binT_irr{i,j} = nanmean(dummy_median_irr, 2);
        binQ25_binT_ctl{i,j}    = nanmean(dummy_Q25_ctl   , 2);
        binQ25_binT_irr{i,j}    = nanmean(dummy_Q25_irr   , 2);
        binQ75_binT_ctl{i,j}    = nanmean(dummy_Q75_ctl   , 2);
        binQ75_binT_irr{i,j}    = nanmean(dummy_Q75_irr   , 2);

    
        % bin all variables according to PRECT
        for k=1:nens
            
            % get PRECT bins (use "tabulate(bin)" for some bin information)
            binP_ctl = ceil(nbins * tiedrank(srex_ctl{7,j,k}) / length(srex_ctl{7,j,k})); % bin according to PRECT (precipitation)
            binP_irr = ceil(nbins * tiedrank(srex_irr{7,j,k}) / length(srex_irr{7,j,k})); % bin according to PRECT (precipitation)
            
            % binning
            [dummy_mean_ctl(:,k), dummy_median_ctl(:,k), ~, dummy_Q25_ctl(:,k), dummy_Q75_ctl(:,k)] = mf_bin(srex_ctl{i,j,k}, binP_ctl, nbins); %#ok<*SAGROW>
            [dummy_mean_irr(:,k), dummy_median_irr(:,k), ~, dummy_Q25_irr(:,k), dummy_Q75_irr(:,k)] = mf_bin(srex_irr{i,j,k}, binP_irr, nbins);
            
        end        
        binmean_binP_ctl{i,j}   = nanmean(dummy_mean_ctl  , 2);
        binmean_binP_irr{i,j}   = nanmean(dummy_mean_irr  , 2);
        binmedian_binP_ctl{i,j} = nanmean(dummy_median_ctl, 2);
        binmedian_binP_irr{i,j} = nanmean(dummy_median_irr, 2);
        binQ25_binP_ctl{i,j}    = nanmean(dummy_Q25_ctl   , 2);
        binQ25_binP_irr{i,j}    = nanmean(dummy_Q25_irr   , 2);
        binQ75_binP_ctl{i,j}    = nanmean(dummy_Q75_ctl   , 2);
        binQ75_binP_irr{i,j}    = nanmean(dummy_Q75_irr   , 2);
    
    
    
    end
end



