%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Driving script for BayesSplit
addpath('Custom_functions_GS')
addpath('CircStat2012a')
addpath('irisFetch')
addpath('FetchData')
addpath('./Reflectivity/toolbox/');
addpath('./Reflectivity/deconvolution_code');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the data to be used
name = 'NoMelt_1layer';

evt_ind      = [ 351,  474,    475,   288,   490   ];%can be array -101.8, -101.4, etc.
target_phase = { 'SKS', 'SKS', 'SKS', 'SKS', 'SKS' };%cell array (inside a cell array
pre          = [ 80,    80,    60,    80,    50    ];%array

%name = 'SynTwoLayers_5evts';
dataName   = 'NoMelt';
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search

Parameters            = make_parameters(name, dataName);
Parameters.parallel   = true;
Parameters.max_layers = 1;

Parameters.r_range             = [ -3.57, 1e-5 ];%mean, std, log
Parameters.f_range             = [ -3.19, 1e-5 ];%mean, std, log

%customize

%for the NoMelt Prior
Parameters.prior_information.use     = false;
Parameters.prior_information.dt      = [ 1 1 ];
Parameters.prior_information.dt_std  = [ 0.33 0.33 ];
Parameters.prior_information.phi     = [ -80 78 ]*pi/180;
Parameters.prior_information.phi_std = [ 15 5 ]*pi/180;
Parameters.prior_information.rot     = [ 5 0 ]*pi/180;
Parameters.prior_information.rot_std = [ 15 5 ]*pi/180;

%%%%%%%%%%%%%%%%%%%%%

BayesSplit
ReportResults
