function [GPSandPosition_table]=GetGPSandPhonePositions(gps_data,depth_data,array_data, parameters, t_serialdate)
% GetGPSandPhonePositions.m is a function that creates a table which
% contains GPS data and hydrophone positions for a given timeframe of the
% encounter (specified in t_serialdate). It assumes that boat moves in a
% straight line (boat positions are smoothed with a 1st order polynomial).
%
%INPUTS
% - gps_data - a table containing GPS information. The following columns are expected:
%       ~ UTC - date and time in UTC;
%       ~ Longitude - boat longitude information;
%       ~ Latitude - boat latitude information.
% - depth_data- a table containing sensor depth information. The following columns are expected:
%       ~ UTC - date and time in UTC;
%       ~ Sensor_0_Depth - sensor 1 depth in m;
%       ~ Sensor_1_Depth - sensor 2 depth in m.
% - array_data - a table containing sensor spacing information.The following columns are expected:
%       ~ Array_name - name of the array that matches the name of the array
%                      specified in parameters struct;
%       ~ distance_behind_boat - array distance behind the boat in m
%       ~ Hyph_1 - distance of the first sensor used from the beginning of the
%                  array in m
% - parameters - a structure containing at least 3 fields:
%               ~ parameters.arrayname - name of the array used to track
%               (should match 'Array_name' in the array_data table)
%               ~ parameters.channels - 1 x 2 vector that indicates which
%               sensors are used for tracking and localization.
%               ~ parameters.d - sensor separation (in m)
% - t_serialdate - a vector of times (in serial date format) for the encounter.
%
%OUTPUTS:
% - GPSandPosition_table - a table containing boat and array position data.
%                           It has the following columns:
%       ~ 'Time_datenum '- date and time in serial date format,
%       ~ 'Time_UTC' - date and time in UTC,
%       ~ 'Boat_Latitude' - boat position - latitude (in decimal degrees),
%       ~ 'Boat_Longitude'- boat position - longitude (in decimal degrees), 
%       ~ 'Boat_pos_x_m'- boat position along x-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%       ~ 'Boat_pos_y_m'- boat position along y-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m), 
%       ~ 'Sensor1_pos_x_m'- sensor 1 position along x-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%       ~ 'Sensor1_pos_y_m'- sensor 1 position along y-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%       ~ 'Sensor1_pos_z_m'- sensor 1 position along z-axis (depth) - Cartesian coordinates 
%                         with respect to start of boat track (in m), 
%       ~ 'Sensor2_pos_x_m'- sensor 2 position along x-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%       ~ 'Sensor2_pos_y_m'- sensor 2 position along y-axis - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%       ~ 'Sensor2_pos_z_m'- sensor 2 position along z-axis (depth) - Cartesian coordinates 
%                         with respect to start of boat track (in m),
%
%
%Pina Gruden, winter 2022, UH Manoa


% datetime array format of times:
Time_UTC=datetime(t_serialdate,'convertfrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');
Time_UTC=Time_UTC(:);

% Get min and max time for the encounter:
t_datetime=datetime(t_serialdate, 'ConvertFrom','datenum');
tmin= t_datetime(1);
tmax= t_datetime(end);

%% Get GPS positions of the BOAT:
% 1) Get indices of times that match tmin&tmax from the overall GPS table:
irow = find(gps_data.UTC(:)>=tmin &  gps_data.UTC(:)<=tmax);
gps_data_select = gps_data([irow(1)-3;irow;irow(end)+3],:); %add more 30 s of data (each point is 10 s of data) at the beginning and end to make sure you cover all possible times

% 2) Interpolate Positions to match time steps used in TDOA tracks:
interp_long= interp1(gps_data_select.UTC,gps_data_select.Longitude,t_datetime);
interp_lat= interp1(gps_data_select.UTC,gps_data_select.Latitude,t_datetime);

% 3) Boat positions in Local Coordinates (X,Y)
start_pos = [interp_lat(1),interp_long(1)]; % start position of the boat
wgs84 = wgs84Ellipsoid;
[x,y,~] = geodetic2enu(interp_lat(:),interp_long(:),0,start_pos(1),start_pos(2),0,wgs84); %Unit is 'meter'
% geodetic2enu - Transform geodetic coordinates to local east-north-up

%check for nan values and remove before fitting a polynomial:
indnan = isnan(x) | isnan(y);
x(indnan) = [];
y(indnan) = [];

% 4) Smooth the boat positions:
p = polyfit(x, y, 1);
y = polyval(p, x);

%% Get positions of the HYDROPHONES:
%//////// DEPTH ////////////
% 1) Get indices of times that match tmin&tmax from the overall DEPTH table:
irow = find(depth_data.UTC(:)>=tmin & depth_data.UTC(:)<=tmax);
depth_data_select = depth_data([irow(1)-1;irow;irow(end)+1],:);

% 2) Interpolate Depths to match time steps used in TDOA tracks:
interp_depth_s1= interp1(depth_data_select.UTC,depth_data_select.Sensor_0_Depth,t_datetime);
interp_depth_s2= interp1(depth_data_select.UTC,depth_data_select.Sensor_1_Depth,t_datetime);

%///////// X,Y Positions (relative to the boat) //////////
row_number=strcmp(array_data.Array_name,parameters.arrayname);
varnames = array_data.Properties.VariableNames;
[~, column_number_a] = ismember('distance_behind_boat', varnames);
match = strfind(varnames, num2str(parameters.channels(1))); %indicates which is the first sensor used
column_number_s1  = ~cellfun('isempty', match);

d_boat_m = array_data{row_number,column_number_a};
 %array distance behind the boat
 %first sensor distance behind boat:
d_s1_m =d_boat_m + array_data{row_number,column_number_s1};
 %second sensor distance behind boat:
d_s2_m = d_s1_m + parameters.d; % parameters.d = sensor separation

%Get relative sensor positions
d_xy = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
%sensor 1 coordinates:
x_s1 = x(2:end) + d_s1_m.*(x(1:end-1)-x(2:end))./d_xy;
x_s1 = [x(1) + d_s1_m.*(x(1)-x(2))/d_xy(1);x_s1];
y_s1 = y(2:end) + d_s1_m.*(y(1:end-1)-y(2:end))./d_xy;
y_s1 = [y(1) + d_s1_m.*(y(1)-y(2))/d_xy(1);y_s1];
%sensor 2 coordinates:
x_s2 = x(2:end) + d_s2_m.*(x(1:end-1)-x(2:end))./d_xy;
x_s2 = [x(1) +  d_s2_m.*(x(1)-x(2))/d_xy(1);x_s2];
y_s2 = y(2:end) + d_s2_m.*(y(1:end-1)-y(2:end))./d_xy;
y_s2 = [y(1) + d_s2_m.*(y(1)-y(2))/d_xy(1);y_s2];


% Create a Table
GPSandPosition_table = table(t_serialdate, Time_UTC,...
    interp_lat,interp_long,x,y,x_s1, y_s1, interp_depth_s1,x_s2, y_s2,interp_depth_s2, ...
    'VariableNames',{'Time_datenum','Time_UTC','Boat_Latitude','Boat_Longitude', ...
    'Boat_pos_x_m','Boat_pos_y_m', 'Sensor1_pos_x_m','Sensor1_pos_y_m','Sensor1_pos_z_m', ...
    'Sensor2_pos_x_m','Sensor2_pos_y_m','Sensor2_pos_z_m'});


end
