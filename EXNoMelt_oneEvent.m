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
name = 'NoMelt_288';

evt_ind      = [ 490   ];%can be array -101.8, -101.4, etc.
target_phase = { 'SKS' };%cell array (inside a cell array
pre          = [ 50    ];%array

%evt_ind      = [ 351,  474,    475,   288,   490,   79    ];%can be array -101.8, -101.4, etc.
%target_phase = { 'SKS', 'SKS', 'SKS', 'SKS', 'SKS', 'PcS' };%cell array (inside a cell array
%pre          = [ 80,    80,    60,    80,    50,    50    ];%array

%name = 'SynTwoLayers_5evts';
dataName   = 'NoMelt';
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%configure the search

Parameters                   = make_parameters(dataName);
Parameters.parallel          = true;
Parameters.get_errors        = false;
Parameters.use_orientations  = false;%include orientations or not in model vector. 

%customize
Parameters.max_layers          = 1;

%for the NoMelt Prior
%%%%%%%%%%%%%%%%%%%%%
BayesSplit
ReportResults
