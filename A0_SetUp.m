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

