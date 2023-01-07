% READ ALL RELEVANT GPS AND SENSOR DATA and SAVE as a TABLE
% This script reads all relevant GPS data for a given encounter from a
% larger GPS database, and also computes sensor positions for each time
% step. Results are saved as .csv table

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had specified paths in:
% a) specify_paths.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

clear, close all
addpath AS_Localization_code/
%//////////////////////////////////////////////////////////////////////////

% --------------------------- Get folders:----------------------------
[folder] = specify_paths;

if ~isempty(folder.rawgps)
    %---------------------LOAD DATA---------------------------
    % LOAD all tables in the raw gps data folder (array, sensor depth and
    % gps info):
    files =dir(fullfile(folder.rawgps,'*.csv'));
    fname   = {files.name};
    fname(~cellfun('isclass', fname, 'char')) = {''};  % Care for non-strings
    
    % Load table containing GPS data
    matchC = reshape(strfind(fname, 'gps'), size(files));
    match  = ~cellfun('isempty', matchC);
    gps_data=readtable([files(match).folder,'/',files(match).name]);

    % Load table containing Depth data
    matchC = reshape(strfind(fname, 'Depth'), size(files));
    match  = ~cellfun('isempty', matchC);
    depth_data=readtable([files(match).folder,'/',files(match).name]);

    % Load table containing Array data:
    matchC = reshape(strfind(fname, 'Array'), size(files));
    match  = ~cellfun('isempty', matchC);
    array_data=readtable([files(match).folder,'/',files(match).name]);

    % LOAD extracted TDOA Tracks, parameters, and time data (from TDOA 
    % tracking master package):
    files =dir(fullfile(folder.tdoas,'*.mat'));
    load([folder.tdoas,files.name],'Tracks','parameters')
    files =dir(fullfile(folder.crosscorr,'*.mat'));
    load([folder.crosscorr,files(1).name],'t_serialdate') % the time vector is the same regardless which file from the same encounter we load
    t_serialdate=t_serialdate(:);

    %-------------------------------------------------------------
    %---------------------CREATE A TABLE---------------------------
    [GPSandPosition_table]=GetGPSandPhonePositions(gps_data,depth_data, ...
        array_data, Tracks, parameters, t_serialdate);

    %-------------------------------------------------------------
    %---------------------SAVE TABLE---------------------------
    writetable(GPSandPosition_table,[folder.gps,'GPSandPosition_table_test.csv']);

    %-------------------------------------------------------------

end