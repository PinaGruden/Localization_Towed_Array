function [GPSandPosition_table]=GetGPSandPhonePositions
% LOAD and Prepare GPS data and corresponding Hydrophone positions

% Load appropriate table containing GPS data
gps_data=readtable('/Users/pinagruden/Dropbox/Pina/HAWAII/Updated_Database_NOAA_FKW/AllLasker_gpsData.csv');

% Load appropriate table containing Depth data
depth_data=readtable('/Users/pinagruden/Dropbox/Pina/HAWAII/Updated_Database_NOAA_FKW/AllLasker_Hydrophone_Depth_Data.csv');

% Load table containing Array data:
array_data=readtable('/Users/pinagruden/Dropbox/Pina/Git_projects/TDOA_tracking_master/Array_Info.csv');

% Load extracted TDOA tracks (to get the times)
load(['/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/' ...
    'Tracking_Towed_Array/LocalizationTest/LaskerAC109_Results/' ...
    'Lasker_AC109_Results.mat'], 'Tracks','parameters')

% Load timevector (from cross-correlograms):
load(['/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/' ...
    'Tracking_Towed_Array/LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/' ...
    'Lasker_AC109_clicks_rawCrossCorrelogram_ALL.mat'],...
    't_serialdate')
t_serialdate=t_serialdate(:);

% datetime array format of times:
Time_UTC=datetime(t_serialdate,'convertfrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');
Time_UTC=Time_UTC(:);

% Get min and max time for the encounter:
tmin = min([Tracks.time_local]); %Tracks.time_local is in datenum format
tmax = max([Tracks.time_local]);

%% Get GPS positions of the BOAT:
% 1) Get indices of times that match tmin&tmax from the overall GPS table:
irow = find(datenum(gps_data.UTC(:))>tmin & datenum(gps_data.UTC(:))<tmax);
gps_data_select = gps_data([irow(1)-1;irow;irow(end)+1],:); %add one more data point at the beginning and end to make sure you cover all possible times

% 2) Interpolate Positions to match time steps used in TDOA tracks:
interp_long= interp1(datenum(gps_data_select.UTC),gps_data_select.Longitude,t_serialdate)';
interp_lat= interp1(datenum(gps_data_select.UTC),gps_data_select.Latitude,t_serialdate)';

% 3) Boat positions in Local Coordinates (X,Y)
start_pos = [interp_lat(1),interp_long(1)]; % start position of the boat
wgs84 = wgs84Ellipsoid;
[x,y,~] = geodetic2enu(interp_lat(:),interp_long(:),0,start_pos(1),start_pos(2),0,wgs84); %Unit is 'meter'


%% Get positions of the HYDROPHONES:
%//////// DEPTH ////////////
% 1) Get indices of times that match tmin&tmax from the overall DEPTH table:
irow = find(datenum(depth_data.UTC(:))>tmin & datenum(depth_data.UTC(:))<tmax);
depth_data_select = depth_data([irow(1)-1;irow;irow(end)+1],:);

% 2) Interpolate Depths to match time steps used in TDOA tracks:
interp_depth_s1= interp1(datenum(depth_data_select.UTC),depth_data_select.Sensor_0_Depth,t_serialdate)';
interp_depth_s2= interp1(datenum(depth_data_select.UTC),depth_data_select.Sensor_1_Depth,t_serialdate)';

%///////// X,Y Positions (relative to the boat) //////////
d_boat_m = array_data{strcmp(array_data.Array_name,parameters.arrayname),4};
 %array distance behind the boat
 %first sensor distance behind boat:
d_s1_m =d_boat_m + array_data{strcmp(array_data.Array_name,parameters.arrayname),5};
 %second sensor distance behind boat:
d_s2_m = d_s1_m + parameters.d; % parameters.d = sensor separation

%Get relative sensor positions
d_xy = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
%sensor 1 coordinates:
x_s1 = x(2:end) + d_s1_m.*(x(1:end-1)-x(2:end))./d_xy;
x_s1 = [0;x_s1];
y_s1 = y(2:end) + d_s1_m.*(y(1:end-1)-y(2:end))./d_xy;
y_s1 = [0;y_s1];
%sensor 2 coordinates:
x_s2 = x(2:end) + d_s2_m.*(x(1:end-1)-x(2:end))./d_xy;
x_s2 = [0;x_s2];
y_s2 = y(2:end) + d_s2_m.*(y(1:end-1)-y(2:end))./d_xy;
y_s2 = [0;y_s2];
% This gets the sensor positions oscillating a bit around the boat
% trajectory. Maybe need to take longer steps than one to another to
% smooth a bit??



% Create a Table
GPSandPosition_table = table(t_serialdate, Time_UTC,...
    interp_lat,interp_long,x,y,x_s1, y_s1, interp_depth_s1,x_s2, y_s2,interp_depth_s2, ...
    'VariableNames',{'Time_datenum','Time_UTC','Boat_Latitude','Boat_Longitude', ...
    'Boat_pos_x_m','Boat_pos_y_m', 'Sensor1_pos_x_m','Sensor1_pos_y_m','Sensor1_pos_z_m', ...
    'Sensor2_pos_x_m','Sensor2_pos_y_m','Sensor2_pos_z_m'});
% save
writetable(GPSandPosition_table,'GPSandPosition_table.csv')

end