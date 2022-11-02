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

%3, 4, 5, 6, 8, 18 (no traces?), 20, (double), 26(double), 28, 29, and 30 bad

evt_ind      = [ 22, 30, 44, 63, 85, 86, 87, 92, 102, 145, 179, 247, 282, ...
    290, 306, 318, 324, 353, 443, 452, 454, 457, 462, 487, 488, 498, 505, ...
    528, 539, 561, 564, 568, 569, 604, 613, 617, 622, 636, 624, 668, 677, ...
    717, 741, 765];

target_phase = { 'SKKS', 'SKS',  'SKKS',  'SKS', 'SKS', 'SKS', 'SKS',  ...
    'SKS', 'SKS', 'SKS', 'SKS', 'SKKS', 'SKS', 'SKS', 'SKKS', 'SKS',   ...
    'SKS', 'SKS', 'SKS', 'SKS', 'SKS', 'SKS', 'SKKS', 'SKS', 'SKKS',   ...
    'SKKS', 'SKS', 'SKS', 'SKS', 'SKS', 'SKS', 'SKKS', 'SKS',  'SKS', 'SKKS', ...
    'SKS', 'SKKS', 'SKS' };

pre          = [ 25,     25,     25,      30,    30,    25,    30,   ...
    30,    15,    10,    25,    25,     15,    15,    25,     10,    ...
    10,    25,    15,    25,    15,    15,    12.5,     15,    15,     ...
    25,     0,    25,    30,    30,    15,     15,     15,     10,  10, ...
    5,    15,     2 ];

station = 1;%[] for all

%name = 'SynTwoLayers_5evts'; 
%dataName   = 'SynTwoLayers_5evts';
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters             = make_parameters('MEEK', 'MEEK');
Parameters.total_time  = 50;%in s

Parameters.get_errors          = false;
Parameters.dtstd               = 0.5;

Parameters.max_layers          = 1;

%Parameters.solver_printout   = true;%show progression of solver or not. Different for different solvers. 

Parameters.sample_rate = 1;%larger by at least 2
Parameters.use_covarience    = false;%include covarience parameters in model vector. 

%%%%%%%%%%%%%%%%%%%%%

BayesSplit
ReportResults
