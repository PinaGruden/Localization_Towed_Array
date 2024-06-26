% SET UP 
% This script reads all relevant data into the Matlab workspace, it also
% selects TDOA tracks and relevant hydrophone positions that will be used
% to obtain final localizations.

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had specified paths and parameters first in:
% a) specify_paths.m
% b) specify_parameters.m

% and if GPS data avaiable, make sure you have run A0_GetGPSData.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

clear, close all
addpath AS_Localization_code/
%//////////////////////////////////////////////////////////////////////////
%% //////////////// 1) Get folders, load data, get parameters /////////////

% --------------------------- a) Get folders:----------------------------
[folder, folder2save2] = specify_paths;
%
% ----------------------- b) Load data from folders: ----------------------
%   i) Load all Pamguard tdoas from detected whistles and clicks (if
%   available)
if ~isempty(folder.pamguard)
    files =dir(fullfile(folder.pamguard,'*.mat'));
    % elimanate OS-X generated files
    idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
    files=files(idx);
    N=size(files,1);
    for k=1:N
        load([folder.pamguard,files(k).name],'All_data_*','parameters')
        pamguard_parameters=parameters;
    end
    
end

%   ii) Load data/tracks extracted by the TDOA tracking package:
%       Load TDOA tracks:
files =dir(fullfile(folder.tdoas,'*.mat'));
% elimanate OS-X generated files
idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
files=files(idx);
for k=1:size(files,1)
    load([folder.tdoas,files(k).name],'Tracks','parameters','scalar_*')
end
% Identify min, max times for these tracks & create a time vector:
timevec=min(vertcat(Tracks(:).time)):parameters.dt:max(vertcat(Tracks(:).time));
parameters.Ntsteps=numel(timevec); %number of time steps

%       Load Cross-correlograms for plotting (if available):
if ~isempty(folder.crosscorr)
    files =dir(fullfile(folder.crosscorr,'*.mat'));
    % elimanate OS-X generated files
    idx=  cellfun(@(x) ~strcmp(x(1),'.'), {files.name}); %get indices of all files that are not '.'
    files=files(idx);
    N=size(files,1);
    Rxy_envelope_both=cell(1,N);
    for k=1:N
        load([folder.crosscorr,files(k).name],'Rxy_envelope_ALL','lags','t_serialdate')
        if ~isempty(strfind(files(k).name,'clicks'))
            Rxy_envelope_both{k}=Rxy_envelope_ALL.*scalar_clicks;
        end
        if ~isempty(strfind(files(k).name,'whistles'))
            Rxy_envelope_both{k}=Rxy_envelope_ALL.*scalar_whistles;
        end
        clear 'Rxy_envelope_ALL'
    end
end

%   iii) Check that the TDOA tracking and Pamguard used the same sensor
%   separation (if applicable)
if ~isempty(folder.pamguard)
    if ~isequal(parameters.d,pamguard_parameters.d)
        error(['The sensor spacing between TDOA tracking and Pamguard detections' ...
            ' does not match. Re-check sensor separation and compute TDOAs for' ...
            ' Pamguard for the sensor separation that matches parameters.d.'])
    end
end

%   iv) Load GPS and array info (If available)
if ~isempty(folder.gps) % Real GPS data is available
    GPSandPosition_table=readtable([folder.gps,'GPSandPosition_table.csv']);
    parameters.get_hyph_pos = 2;
    parameters.GPSTable=GPSandPosition_table;
else
    parameters.get_hyph_pos = 1; % Real GPS data is not available
end

% ---------------------- c) Get paramters: ----------------------
[AS_params,BA_params] = specify_parameters(parameters);


%//////////////////////////////////////////////////////////////////////////
%% ///////// 2) Get All Relevant TDOA tracks around the beam //////////////

% these tracks are selected based on angles around the beam that you 
% specified in the specify_parameters.m 

%-------- a) Search for all tracks within x degrees of the beam -----------
indx=nan(size(Tracks,2),1);
for k=1:size(Tracks,2)
indx(k)=ismembertol(0,Tracks(k).tdoa,AS_params.tdoa_cutoff); 
end
indx=logical(indx);
Nsources=sum(indx);
Tracks_selected=Tracks(indx);
clear indx

%-------------- c) PLOT All tracks and Selected tracks:------------------

% Plot against cross-correlograms and against Pamguard detections (if
% available)
if ~isempty(folder.pamguard) % If Pamguard detections are available
    switch parameters.signal_type
        case 'both'
            plot_tracks(folder,Tracks,Tracks_selected,t_serialdate,lags,parameters, ...
                AS_params.tdoa_cutoff,Rxy_envelope_both,All_data_c,All_data_w)
        case 'whistles'
            plot_tracks(folder,Tracks,Tracks_selected,t_serialdate,lags,parameters, ...
                AS_params.tdoa_cutoff,Rxy_envelope_both,All_data_w)
        case 'clicks'
            plot_tracks(folder,Tracks,Tracks_selected,t_serialdate,lags,parameters, ...
                AS_params.tdoa_cutoff,Rxy_envelope_both,All_data_c)
    end
else % If Pamguard detections are not available
    plot_tracks(folder,Tracks,Tracks_selected,t_serialdate,lags,parameters, ...
        AS_params.tdoa_cutoff,Rxy_envelope_both)
end





