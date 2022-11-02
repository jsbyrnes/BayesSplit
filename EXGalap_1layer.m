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

evt_ind      = [ 1232,   1305,   1450,    639,   979,   1208,  1500, ...
    1505,  687,   1620,  1712,  1599,   1586,  1534,  1532,   1520,  ...
    1370,  1227,  1642,  687,   1512,  1101,  1101,   1100,  1100,   ...
    1100,   1038,  1635,  1208,  756,   1412,   1219,   1219,   741, ...
    741, 605,   605,    1144 ];

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

evt_ind     ([ 3, 4, 5, 6, 8, 18, 20, 26, 28, 29, 30 ]) = [];
target_phase([ 3, 4, 5, 6, 8, 18, 20, 26, 28, 29, 30 ]) = [];
pre         ([ 3, 4, 5, 6, 8, 18, 20, 26, 28, 29, 30 ]) = [];

%name = 'SynTwoLayers_5evts';
%dataName   = 'SynTwoLayers_5evts';
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters             = make_parameters('IGUANA-1layer', 'IGUANA');
Parameters.parallel    = false;
Parameters.total_time  = 50;%in s
Parameters.max_layers  = 1;

Parameters.dtstd               = 0.5;
Parameters.use_tSA           = true;

%Parameters.sample_rate = 1;%larger by at least 2

%%%%%%%%%%%%%%%%%%%%%

BayesSplit
ReportResults
