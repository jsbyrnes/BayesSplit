clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Driving script for BayesSplit
addpath('Custom_functions_RF')
addpath('CircStat2012a')
addpath('irisFetch')
addpath('FetchData')
addpath('./Reflectivity/toolbox/')
addpath('./Reflectivity/deconvolution_code')
%addpath('./fastnonlinearACD/')
%addpath('./ParallelFastNonLinearACD-master/')
%addpath('./fminlbfgs_version2c')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the data to be used 
name     = 'HLP_72RF';
dataName = 'HLP_small4';

evt_ind      = 71;
target_phase = { 'P' };
pre          = 10;
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters                 = make_parameters(name, 'HLP_small4');
Parameters.parallel        = true;
Parameters.get_errors      = false;
Parameters.total_time      = 50;%in s
Parameters.solver_printout = true;%show progression of solver or not. Different for different solvers. 
% Parameters.low_pass        = 0.25;
% Parameters.sample_rate     = 3*0.25;%larger by at least 2

Parameters.reset_rounds      = 10;
Parameters.max_layers        = 1;
Parameters.rotation_std      = 1e-2*pi/180;

Parameters.low_pass    = 2.5;
Parameters.sample_rate = 5;%larger by at least 2

Parameters.high_pass   = 1/5;%3/(Parameters.total_time);%in seconds

Parameters.vertical_names    = { 'BHZ', 'HHZ' };
Parameters.max_gaussians     = 20; 
Parameters.use_orientations  = false;%include orientations or not in model vector. 
Parameters.use_polarization  = false;
Parameters.use_covarience    = true;
Parameters.wavelet           = true;

Parameters.r_range             = [ log(0.1), 0.5 ];%mean, std, log
Parameters.f_range             = [ log(1), 0.5 ];%mean, std, log
%%%%%%%%%%%%%%%%%%%%%

BayesSplit