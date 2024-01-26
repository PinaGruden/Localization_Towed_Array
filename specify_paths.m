function [folder, folder2save2] = specify_paths
% specify_paths.m specifies paths for "Localization_Towed_Array" package
%
%CONFIGURE PATHS to :
% 1) Results of Extracted TDOA tracks (from TDOA_tracking_masters package)
% 2) Cross-correlogram info (from TDOA_tracking_masters package)
% 3) Results of Extracted Pamguard detections of clicks and whistles (if 
% available,from Extract_Pamguard_detections package)
% 4) Folder where all raw GPS and Array data is
% 5) Folder where GPS and Sensor Position table should be saved
% 6) Folder where localization results should be saved to
%
% OUTPUT:
% - folder - a structure specifying paths to where data is
% - folder2save2 - a structure specifying paths to where results are saved to.
%
%Pina Gruden, 2023, UH Manoa


%1) Path to Results of Extracted TDOA tracks (TDOA_tracking_masters package):
myfolder = '/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Troubleshooting_Tracking_and_Localization_Packages/Issue_33/3_Localize/Data/TDOA_tracks/';
if not(isfolder(myfolder)) % If the folder doesnt exist - throw error since you need data to localize.
error(['The folder ', myfolder, ' does not exists. There is no data ' ...
    'for localization! Run "TDOA_tracking_master" package first!'])
end
s=what(myfolder); 
folder.tdoas = [s.path,'/'];

%2) Path to Cross-correlogram information (TDOA_tracking_masters package):
myfolder = './Test_example/Data/Crosscorrelograms/';
if not(isfolder(myfolder)) % If the folder doesnt exist
    mkdir(myfolder) %make a folder    
    disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
end
s=what(myfolder); 
folder.crosscorr = [s.path,'/'];

%3) Path to Results of Extracted Pamguard detections (Extract_Pamguard_detections package):
% If they are not available specify: myfolder = [];
myfolder='./Test_example/Data/Pamguard_detections/';
if ~isempty(myfolder)
    if not(isfolder(myfolder)) % If the folder doesnt exist
        mkdir(myfolder) %make a folder
        disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
    end
    s=what(myfolder);
    folder.pamguard = [s.path,'/'];
else
    folder.pamguard = [];
end

% 4) Path to all raw GPS and Array data
% If GPS data not available specify: myfolder = [];
myfolder='./Test_example/Data/Raw_GPSandArray_info/';
if ~isempty(myfolder)
    if not(isfolder(myfolder)) % If the folder doesnt exist
        mkdir(myfolder) %make a folder
        disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
    end
    s=what(myfolder);
    folder.rawgps = [s.path,'/'];
else
    folder.rawgps = [];
end

% 5) Path to GPS and Sensor Position table will be stored:
% If GPS data not available specify: myfolder = [];
myfolder='./Test_example/Data/GPSandPos_info/';
if ~isempty(myfolder)
    if not(isfolder(myfolder)) % If the folder doesnt exist
        mkdir(myfolder) %make a folder
        disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
    end
    s=what(myfolder);
    folder.gps = [s.path,'/'];
else
    folder.gps = [];
end

% 6) Path to where Localization Results will be stored:
myfolder='/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Troubleshooting_Tracking_and_Localization_Packages/Issue_33/3_Localize/Results/';
if not(isfolder(myfolder)) % If the folder doesnt exist
    mkdir(myfolder) %make a folder    
    disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
end
s=what(myfolder); 
folder2save2 = [s.path,'/'];


end