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
    % elimanate OS-X generated files
    idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
    files=files(idx);
    fname   = {files.name};
    fname(~cellfun('isclass', fname, 'char')) = {''};  % Care for non-strings
    
    % Load table containing GPS data
    matchC = reshape(strfind(lower(fname), 'gps'), size(files));
    match  = ~cellfun('isempty', matchC);
    gps_data=readtable([files(match).folder,'/',files(match).name]);

    % Load table containing Depth data
    matchC = reshape(strfind(lower(fname), 'depth'), size(files));
    match  = ~cellfun('isempty', matchC);
    depth_data=readtable([files(match).folder,'/',files(match).name]);

    % Load table containing Array data:
    matchC = reshape(strfind(lower(fname), 'array'), size(files));
    match  = ~cellfun('isempty', matchC);
    array_data=readtable([files(match).folder,'/',files(match).name]);

    % LOAD extracted TDOA Tracks, parameters, and time data (from TDOA 
    % tracking master package):
    files =dir(fullfile(folder.tdoas,'*.mat'));
    % elimanate OS-X generated files
    idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
    files=files(idx);
    load([folder.tdoas,files.name],'parameters')
    files =dir(fullfile(folder.crosscorr,'*.mat'));
    % elimanate OS-X generated files
    idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
    files=files(idx);
    load([folder.crosscorr,files(1).name],'t_serialdate') % the time vector is the same regardless which file from the same encounter we load
    t_serialdate=t_serialdate(:);

    %-------------------------------------------------------------
    %---------------------CREATE A TABLE---------------------------
    [GPSandPosition_table]=GetGPSandPhonePositions(gps_data,depth_data, ...
        array_data,parameters, t_serialdate);

    %-------------------------------------------------------------
    %---------------------SAVE TABLE---------------------------
    writetable(GPSandPosition_table,[folder.gps,'GPSandPosition_table.csv']);

    %-------------------------------------------------------------

end