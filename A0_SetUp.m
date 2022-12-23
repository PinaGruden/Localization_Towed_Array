% SET UP 

% ///////////////// 1) Get folders, load data, get parameters /////////////
% CHANGE BELOW - copy over the function to this project & modify 
% Add path for plotting- change to reflect your path to the TDOA tracking package:
addpath('/Users/pinagruden/Dropbox/Pina/Git_projects/TDOA_tracking_master'); 

% --------------------------- a) Get folders:----------------------------
[folder, folder2save2] = specify_paths;

% ----------------------- b) Load data from folders: ----------------------
%   i) Load all Pamguard tdoas from detected whistles and clicks (if
%   available)
if ~isempty(folder.pamguard)
    files =dir(fullfile(folder.pamguard,'*.mat'));
    for k=1:size(files,1)
        load([folder.pamguard,files(k).name],'All_data_*','parameters')
        pamguard_parameters=parameters;
    end
end
%   ii) Load data/tracks extracted by the TDOA tracking package:
load([folder.tdoas,'*.mat'],'Tracks','parameters','scalar_clicks','scalar_whistles')
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/Tracking_Towed_Array/LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_clicks_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL')
Rxy_envelope_both{1}=Rxy_envelope_ALL.*scalar_clicks;
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/Tracking_Towed_Array/LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_whistles_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL','lags','t_serialdate')
Rxy_envelope_both{2}=Rxy_envelope_ALL.*scalar_whistles;

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
GPSandPosition_table=readtable('GPSandPosition_table.csv');

% ---------------------- c) Get paramters: ----------------------
[AS_params,BA_params] = specify_parameters(parameters);

