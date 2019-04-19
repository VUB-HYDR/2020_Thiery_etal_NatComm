

% --------------------------------------------------------------------
% function to compute LRs and related diagnostics
% --------------------------------------------------------------------


function [prob_ref2exp, LR_ref2exp, LR_ref2exp_lp, LR_ref2exp_ip, LR_ref2exp_pinb, ...
          boxplotdata_LR_ref2exp_lp, boxplotdata_LR_ref2exp_ip, boxplotdata_LR_ref2exp_pinb, LR_ref2exp_srex_lp, LR_ref2exp_srex_ip] = mf_LR(perc_ref, perc_ref2exp, island, isirr, ispinb, issrex)



% --------------------------------------------------------------------
% Manipulations
% --------------------------------------------------------------------


% get number of probabilities(=percentages)
npercentages = size(perc_ref,3);


% get number of srex regions
nreg = size(issrex,3);


% convert from percentage to probability
prob_ref     = 1 - (perc_ref     ./ 100);
prob_ref2exp = 1 - (perc_ref2exp ./ 100);


% compute risk ratio (LR) maps
LR_ref2exp = prob_ref2exp ./ prob_ref;


% get spatial median LRs
[~, LR_ref2exp_lp, LR_ref2exp_ip, LR_ref2exp_pinb] = mf_fieldmedian(LR_ref2exp, island, isirr, ispinb);


% prepare boxplotdata
for i=1:npercentages
    
    % isolate percentage map
    LR_ref2exp_i                     = LR_ref2exp(:,:,i);
    
    % generate boxplotdata - all land
    boxplotdata_LR_ref2exp_lp(:,i) = LR_ref2exp_i(island); %#ok<*AGROW>
    
    % generate boxplotdata - irrigated land
    boxplotdata_LR_ref2exp_ip(:,i) = LR_ref2exp_i(isirr);

    % generate boxplotdata - PINB
    boxplotdata_LR_ref2exp_pinb(:,i) = LR_ref2exp_i(ispinb);
    
end


% compute median LRs and FARs per srex region
for i=1:nreg
    
    % get LRs
    [~, LR_ref2exp_srex_lp(:,i), LR_ref2exp_srex_ip(:,i)] = mf_fieldmedian(LR_ref2exp, issrex(:,:,i) & island, issrex(:,:,i) & isirr); 

end


% compute Fraction of Attributable Risk (FAR) maps
%FAR_ref2exp = 1 - ( prob_ref ./ prob_ref2exp );


% get spatial median FARs
%[~, FAR_ref2exp_lp, FAR_ref2exp_ip, FAR_ref2exp_pinb] = mf_fieldmedian(FAR_ref2exp, island, isirr, ispinb);


end

