
% --------------------------------------------------------------------
% main script to postprocess and visualise CESM output: CTL and IRR
% For 'event D&A' study
% --------------------------------------------------------------------



tic


% clean up
clc;
clear;
close all;


% flags
flags.surv  = 0; % 0: do not compute survival plots
                 % 1: compute survival plots
flags.LR    = 0; % 0: do not compute likelihood ratio (LR)
                 % 1: compute likelihood ratio (LR)
flags.perc  = 0; % 0: do not compute percentile changes
                 % 1: compute percentile changes
flags.valp  = 0; % 0: do not compute values used in the paper
                 % 1: compute values used in the paper
flags.refr  = 0; % 0: reference run includes irrigation (20ci) - default
                 % 1: reference run excludes irrigation (20cc) - only for sfig5
flags.plot  = 1; % 0: do not plot
                 % 1: plot



% --------------------------------------------------------------------
% initialisation
% --------------------------------------------------------------------


% declare globals
global island                                                              %#ok<NUSED>


% add matlab scripts directory to path
addpath(genpath('C:\Users\u0079068\Documents\Research\matlab_scripts'));


% add directory containing nc files to path
addpath(genpath('C:\Users\u0079068\Documents\Research\CESM_present\ncfiles'));
addpath(genpath('C:\Users\u0079068\Documents\Research\CESM_present_DA\ncfiles'));


% initialise model parameters
nens        = 5;   % number of ensemble members


% initialise time parameters - CESM
time_begin  = [1981, 1, 1, 0,0,0];
time_end    = [2010,12,31,23,0,0];
years       = (time_begin(1):time_end(1))';
nyears      = length(years);


% initialise srex variable/unit/region names                           
srex_vars    = {'QIRRIG'  , 'TSA'  , 'Qle'    , 'FSH'    , 'FSA'     , 'FIRA'    , 'PRECT'        , 'QRUNOFF'      , 'TREFHT', 'TREFHTMN'  , 'TREFHTMX'  }; % srex variables considered in this study ...
srex_ylabels = {'Q_i_r_r' , 'T_2_m', 'LHF'    , 'SHF'    , 'SW_n_e_t', 'LW_n_e_t', 'Precipitation', 'Q_R_U_N_O_F_F', 'T_2_m' , 'T_{2m,min}', 'T_{2m,max}'}; % ... their y-axis label,
srex_units   = {'mm d^-^1', 'K'    , 'W m^-^1', 'W m^-^1', 'W m^-^1' , 'W m^-^1' , 'mm d^-^1'     , 'mm d^-^1'     , 'K'     , 'K'         , 'K'         }; % ... their units,
srex_block   = {'lnd'     , 'lnd'  , 'lnd'    , 'lnd'    , 'lnd'     , 'lnd'     , 'atm'          , 'lnd'          , 'atm'   , 'atm'       , 'atm'       }; % ... and whether it's a land or atmospheric variable
srex_reg     = {'WNA', 'CNA', 'MED', 'WAS', 'SAS', 'EAS'};                                                                                           % srex regions considered in this study


% initialise percentages (probabilities = 1-percentage/100) used for percentile calculations
percentages = [50.000   90.000   95.000   97.500   99.000   99.500   99.750   99.900];%   99.950   99.975   99.990];



% --------------------------------------------------------------------
% load data
% --------------------------------------------------------------------


ms_load



% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


ms_manip



% --------------------------------------------------------------------
% Model evaluation
% --------------------------------------------------------------------


if flags.surv == 1   
   ms_surv;
end



% --------------------------------------------------------------------
% Statistical significance
% --------------------------------------------------------------------


if flags.LR == 1   
   ms_LR;
end



% --------------------------------------------------------------------
% Percentile change
% --------------------------------------------------------------------


if flags.perc == 1
   ms_perc
end



% --------------------------------------------------------------------
% get values used in the paper
% --------------------------------------------------------------------


if flags.valp == 1
   ms_valp
end



% --------------------------------------------------------------------
% visualise output
% --------------------------------------------------------------------


if flags.plot == 1
   ms_plotscript
end




toc
