
% --------------------------------------------------------------------
% subroutine to compute likelihood ratio
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


% convert from percentage to probability
prob_ref     = 1 - (perc_ref ./ 100);


% call function computing likelihood ratios and related diagnostics - TREFHTMX
[prob_pic2ctl, LR_pic2ctl, LR_pic2ctl_lp, LR_pic2ctl_ip, LR_pic2ctl_pinb, boxplotdata_LR_pic2ctl_lp, boxplotdata_LR_pic2ctl_ip, boxplotdata_LR_pic2ctl_pinb, LR_pic2ctl_srex_lp, LR_pic2ctl_srex_ip] = mf_LR(perc_ref, perc_pic2ctl, island, isirr, ispinb, issrex);
[prob_ctl2irr, LR_ctl2irr, LR_ctl2irr_lp, LR_ctl2irr_ip, LR_ctl2irr_pinb, boxplotdata_LR_ctl2irr_lp, boxplotdata_LR_ctl2irr_ip, boxplotdata_LR_ctl2irr_pinb, LR_ctl2irr_srex_lp, LR_ctl2irr_srex_ip] = mf_LR(perc_ref, perc_ctl2irr, island, isirr, ispinb, issrex);
[prob_pic2irr, LR_pic2irr, LR_pic2irr_lp, LR_pic2irr_ip, LR_pic2irr_pinb, boxplotdata_LR_pic2irr_lp, boxplotdata_LR_pic2irr_ip, boxplotdata_LR_pic2irr_pinb, LR_pic2irr_srex_lp, LR_pic2irr_srex_ip] = mf_LR(perc_ref, perc_pic2irr, island, isirr, ispinb, issrex);
