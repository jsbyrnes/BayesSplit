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
name     = 'HLP_142_cIR';
dataName = 'HLP_small4';

evt_ind      = 142;
target_phase = { 'SKS' };
pre          = 20;
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters                 = make_parameters('HLP_142_tSA', 'HLP_small4');
Parameters.parallel        = false;
Parameters.get_errors      = true;
Parameters.total_time      = 60;%in s
Parameters.solver_printout = true;%show progression of solver or not. Different for different solvers. 
% Parameters.low_pass        = 0.25;
% Parameters.sample_rate     = 3*0.25;%larger by at least 2

Parameters.reset_rounds      = 10;
Parameters.max_layers        = 1;
Parameters.rotation_std      = 1e-2*pi/180;

Parameters.low_pass    = 1/5;
Parameters.sample_rate = Parameters.low_pass*3;%larger by at least 2

Parameters.high_pass   = 1/25;%3/(Parameters.total_time);%in seconds

Parameters.use_orientations  = false;%include orientations or not in model vector. 
Parameters.use_polarization  = false;
Parameters.use_covarience    = true;
Parameters.wavelet           = true;
%%%%%%%%%%%%%%%%%%%%%

BayesSplit
ReportResults