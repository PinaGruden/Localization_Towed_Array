% SET UP 
% This script reads all relevant data into the Matlab workspace, it also
% selects TDOA tracks and relevant hydrophone positions that will be used
% to obtain final localizations.

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had specified paths and parameters first in:
% a) specify_paths.m
% b) specify_parameters.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

clear, close all
addpath AS_Localization_code/
%//////////////////////////////////////////////////////////////////////////
%% //////////////// 1) Get folders, load data, get parameters /////////////

% --------------------------- a) Get folders:----------------------------
[folder, folder2save2] = specify_paths;

% ----------------------- b) Load data from folders: ----------------------
%   i) Load all Pamguard tdoas from detected whistles and clicks (if
%   available)
if ~isempty(folder.pamguard)
    files =dir(fullfile(folder.pamguard,'*.mat'));
    N=size(files,1);
    for k=1:N
        load([folder.pamguard,files(k).name],'All_data_*','parameters')
        pamguard_parameters=parameters;
    end
    
end

%   ii) Load data/tracks extracted by the TDOA tracking package:
%       Load TDOA tracks:
files =dir(fullfile(folder.tdoas,'*.mat'));
for k=1:size(files,1)
    load([folder.tdoas,files(k).name],'Tracks','parameters','scalar_*')
end

%       Load Cross-correlograms for plotting (if available):
if ~isempty(folder.crosscorr)
    files =dir(fullfile(folder.crosscorr,'*.mat'));
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
if ~isempty(folder.gps)
    GPSandPosition_table=readtable([folder.gps,'GPSandPosition_table.csv']);
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

%--- b) Identify min, max times for these tracks & create a time vector---
timevec=min(vertcat(Tracks_selected(:).time)):parameters.dt:max(vertcat(Tracks_selected(:).time));
Ntsteps=numel(timevec); %number of time steps


%-------------- c) PLOT All tracks and Selected tracks:------------------

% Plot against cross-correlograms and against Pamguard detections (if
% available)

plot_tracks(folder,Tracks,Tracks_selected,t_serialdate,lags,parameters, ...
    Rxy_envelope_both,All_data_c,All_data_w)



%//////////////////////////////////////////////////////////////////////////
%% ////////////// 3) Get Relevant Hydrophone positions ///////////////////

switch BA_params.get_hyph_pos
    case 1 % SIMUALTED GPS DATA
        [hyph_pos,x_boat,y_boat] = simulate_array_pos(Ntsteps,BA_params);

    case 2 % REAL GPS DATA
        % Get Hydrophone positions from data
         hyph_pos=nan(2,2,Ntsteps);
        for t=1:Ntsteps
            hyph_pos(:,:,t)= [GPSandPosition_table.Sensor1_pos_x_m(t),...
                GPSandPosition_table.Sensor1_pos_y_m(t);...
                GPSandPosition_table.Sensor2_pos_x_m(t),...
                GPSandPosition_table.Sensor2_pos_y_m(t)]; %[x1,y1; x2,y2];
        end

end






