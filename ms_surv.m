
% --------------------------------------------------------------------
% subroutine to compute survival functions
% note: preferably run "main"
% --------------------------------------------------------------------


% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


% loop over srex variables
for i=1:length(srex_vars) 

    % loop over srex regions
    for j=1:nreg
% for i=11
%     for j=5

        
        % generate 1 time series for the nens ensemble members        
        dummy_ctl = [];
        dummy_irr = [];
        for k=1:nens
            
            % dummy variable
            dummy_ctl = [dummy_ctl srex_ctl{i,j,k}]; %#ok<*AGROW>
            dummy_irr = [dummy_irr srex_irr{i,j,k}];
            
            
        end        
    
    
        % Kaplan-Meier estimate of the cumulative distribution function (cdf), also known as the empirical cdf
        % Plot the survivor function for the data with 99% confidence bounds.
        [surv_f_pic{i,j}, surv_x_pic{i,j}, ~, ~] = ecdf(srex_pic{i,j}, 'function', 'survivor', 'alpha', 0.01, 'bounds', 'on');
        [surv_f_ctl{i,j}, surv_x_ctl{i,j}, ~, ~] = ecdf(dummy_ctl(:) , 'function', 'survivor', 'alpha', 0.01, 'bounds', 'on'); %#ok<*SAGROW>
        [surv_f_irr{i,j}, surv_x_irr{i,j}, ~, ~] = ecdf(dummy_irr(:) , 'function', 'survivor', 'alpha', 0.01, 'bounds', 'on');


        % Kernel smoothing density estimate
        [ksdens_f_pic{i,j}, ksdens_x_pic{i,j}] = ksdensity(srex_pic{i,j});
        [ksdens_f_ctl{i,j}, ksdens_x_ctl{i,j}] = ksdensity(dummy_ctl(:));
        [ksdens_f_irr{i,j}, ksdens_x_irr{i,j}] = ksdensity(dummy_irr(:));
        
        
        % Risk Ratio (LR)
        [LR_ts_pic2ctl{i,j}, LR_ts_ci_pic2ctl{i,j}, FAR_ts_pic2ctl{i,j}] = mf_riskratio(srex_pic{i,j}, dummy_ctl(:), 0.05, percentages);  % change in LR due to GHG emissions
        [LR_ts_ctl2irr{i,j}, LR_ts_ci_ctl2irr{i,j}, FAR_ts_ctl2irr{i,j}] = mf_riskratio(dummy_ctl(:) , dummy_irr(:), 0.05, percentages);  % change in LR due to GHG emissions
        [LR_ts_pic2irr{i,j}, LR_ts_ci_pic2irr{i,j}, FAR_ts_pic2irr{i,j}] = mf_riskratio(srex_pic{i,j}, dummy_irr(:), 0.05, percentages);  % change in LR due to GHG emissions
        
    end
end
