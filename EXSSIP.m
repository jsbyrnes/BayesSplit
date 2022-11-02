clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Driving script for BayesSplit
addpath('Custom_functions_GS')
addpath('CircStat2012a')
addpath('irisFetch')
addpath('FetchData')
addpath('./Reflectivity/toolbox/')
addpath('./Reflectivity/deconvolution_code')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the data to be used 
name     = 'SSIP';
dataName = 'SSIP';

load('SSIP_picks_small')

%evt_ind      = evt_ind([ 17 ]);%[ 1 3 17 19 20 ]
%pre          = pre([ 17 ]);
%target_phase = target_phase([ 17 ]);
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters                 = make_parameters('SSIP', 'SSIP');
Parameters.parallel        = true;
Parameters.get_errors      = false;
Parameters.total_time      = 50;%in s
% Parameters.low_pass        = 0.25;
% Parameters.sample_rate     = 3*0.25;%larger by at least 2

Parameters.reset_rounds      = 2;
Parameters.max_layers        = 1;

Parameters.low_pass    = 0.15;
Parameters.sample_rate = 0.45;%larger by at least 2

Parameters.high_pass   = 1/50;%3/(Parameters.total_time);%in seconds

Parameters.solver_printout   = true;%show progression of solver or not. Different for different solvers. 
Parameters.use_orientations  = true;%include orientations or not in model vector. 
Parameters.use_polarization  = true;
Parameters.use_tSA           = true;
Parameters.use_covarience    = true;
Parameters.wavelet           = true;
%%%%%%%%%%%%%%%%%%%%%

BayesSplit
ReportResults