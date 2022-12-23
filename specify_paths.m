function [folder, folder2save2] = specify_paths
%CONFIGURE PATHS to :
% 1) Results of Extracted TDOA tracks (from TDOA_tracking_masters package)
% 2) Results of Extracted Pamguard detections of clicks and whistles (if 
% available,from Extract_Pamguard_detections package)
% 3) GPS and Sensor Position table
% 4) Folder where localization results should be saved to

%1) Path to Results of Extracted TDOA tracks (TDOA_tracking_masters package):
s=what('./Test_example/Data/TDOA_tracks/'); 
folder.tdoas = [s.path,'/'];

%2) Path to Results of Extracted Pamguard detections (Extract_Pamguard_detections package):
% If they are not avaiable specify folder.pagmuard = [];
s=what('./Test_example/Data/Pamguard_detections/'); 
folder.pagmuard = [s.path,'/'];

% 3) Path to GPS and Sensor Position table (Table obtained by running
% GetGPSandPhonePositions.m)
% If they are not avaiable specify folder.gps = [];
s=what('./Test_example/Data/GPSandPos_info/'); 
folder.gps = [s.path,'/'];


% 4) Path to where Results will be stored:
s=what('./Test_example/Results/'); 
folder2save2 = [s.path,'/'];


end