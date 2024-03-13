function [AS_params,BA_params] = specify_parameters(parameters)
% Specify parameters for towed array localization:
% Make sure you change/adjust parameters under ALL captions called "Changable"
%
%INPUTS:
% - parameters : a structure containing array and encounter info (used to
%   obtain cross-correlograms)- the minimum fields required are:
%   ~ c : speed of sound - a scalar 
%   ~ d : distance between the two sensors - a scalar
%   ~ dt : how much time elapses in each time step in s - a scalar
%   ~ get_hyph_pos : a scalar indicating whether simulated or real
%                  GPS data is used (1 = simulated; 2 = real)
%       - When get_hyph_pos==2, optional input is also
%           ~ GPSTable: Table containing GPS data 
%   ~ Ntsteps : a scalar indicating how many time steps are in the data
%
%OUTPUTS:
% - AS_params= a structure containing parameters for ambiguity surface 
%              computation. It has the following fields:
%   ~ sig : standard deviation (std) for the Gaussian for ambiguity surface
%   ~ sig_hyperbolas : std for plotting intersecting hyperbolas
%   ~ xrange : x range in m (for grid computation) - a 1x2 vector
%   ~ yrange : y range in m (for grid computation) - a 1x2 vector
%   ~ dx : grid step size in m (resolution) in x direction - a scalar
%   ~ dy : grid step size in m (resolution) in y direction - a scalar
%   ~ c : speed of sound - a scalar
%   ~ mean_horiz_swimspeed : mean horizontal swim speed [m/s] - a scalar
%   ~ bearing_cutoff : cutoff for min bearing around the beam to be
%                      considered in the analysis - a scalar
%   ~ tdoa_cutoff : cutoff for min tdoa around the beam to be
%                      considered in the analysis - a scalar
%   ~ Ngp_x : number of grid points in x - a scalar
%   ~ Ngp_y : number of grid points in y - a scalar
%   ~ X : a grid of X positions - a matrix of dimension Ngp_y x Ngp_x
%   ~ Y : a grid of Y positions - a matrix of dimension Ngp_y x Ngp_x
%   ~ N : number of grid points to evaluate - a scalar
%   ~ wpos : a grid of N x [x,y] coordinates
%
% - BA_params= a structure containing parameters for the boat and array. It
%               has the following fields:
%   ~ get_hyph_pos : a scalar indicating whether simulated or real
%                  GPS data will be used (1 = simulated; 2 = real)
%        - When get_hyph_pos==1 the following fields are included:
%           * mb : line gradient- a scalar- specifies line on which boat
%                  moves 
%           * b : intersect with y axis- a scalar- specifies line on which 
%                 boat moves 
%           * boatspeed: boat speed in m/s
%           * timestep :how much time elapses in each time step in s 
%           * d : distance between the sensors
%           * hyph_pos - sensor positions per time step- a 2 x N x M array, 
%           where N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps.
%           * boat_pos - boat position (in relative coordinates)- M x 2 matrix, 
%           where M is number of time steps, and each row is [x,y] coordinate 
%           for that time step
%           * boat_start_latlong - boat start position in latitude and longitude 
%              (decimal degrees), 1 x 2 vector [latitude, longitude];
%       - When get_hyph_pos==2 the following fields are included:
%           * timestep :how much time elapses in each time step in s
%           * hyph_pos - sensor positions per time step- a 2 x N x M array, 
%           where N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps.
%           * boat_pos - boat position (in relative coordinates)- M x 2 matrix, 
%           where M is number of time steps, and each row is [x,y] coordinate 
%           for that time step
%           * boat_start_latlong - boat start position in latitude and longitude 
%              (decimal degrees), 1 x 2 vector [latitude, longitude];
%
%
%Pina Gruden, Dec 2022, UH Manoa

%% //////////////////// SET PARAMETERS //////////////////// 


%///////////////////// PARAMETERS for the BOAT and ARRAY /////////////////

% Determine wether you had real GPS data available or not (this is
% encapsulated in the parameter parameters.get_hyph_pos (obtained in script A0_SetUp.m)
% and is based on wether you had any data in the folder or not).
%  parameters.get_hyph_pos will be either 1 or 2- i.e.
% ~ parameters.get_hyph_pos = 1 will simulate boat movement and hydrophone positions.
% ~ parameters.get_hyph_pos = 2 will use the real data

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~

%ONLY CHANGE BELOW IF YOU ARE SIMULATING GPS POSITIONS - and ONLY CHANGE
%the parameters **INSIDE** the IF STATMENT:
if parameters.get_hyph_pos==1
    boatspeed_kts = 10; % boat speed in knots
    %specify line gradient and intersect along which the boat moves:
    BA_params.mb=-0.5; % line gradient
    BA_params.b=10; %intersect with y axis
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BA_params.get_hyph_pos = parameters.get_hyph_pos; 

switch BA_params.get_hyph_pos

    case 1 % SIMUALTED GPS DATA
        BA_params.boatspeed= boatspeed_kts/1.944; %boat speed in m/s
        BA_params.timestep = parameters.dt; % how much time elapses in each time step in s 
        BA_params.d = parameters.d; % distance between the sensors
        [hyph_pos,x_boat,y_boat] = simulate_array_pos(parameters.Ntsteps,BA_params);
        BA_params.hyph_pos=hyph_pos;
        boat_pos(:,1)=x_boat;
        boat_pos(:,2)=y_boat;
        BA_params.boat_pos=boat_pos;
        BA_params.boat_start_latlong =[20.4606,-156.9508];

    case 2 % REAL GPS DATA
        GPSandPosition_table=parameters.GPSTable;
        BA_params.timestep = parameters.dt; % how much time elapses in each time step in s
        % Get Hydrophone positions from data
        hyph_pos=nan(2,2,parameters.Ntsteps);
        for t=1:parameters.Ntsteps
            hyph_pos(:,:,t)= [GPSandPosition_table.Sensor1_pos_x_m(t),...
                GPSandPosition_table.Sensor1_pos_y_m(t);...
                GPSandPosition_table.Sensor2_pos_x_m(t),...
                GPSandPosition_table.Sensor2_pos_y_m(t)]; %[x1,y1; x2,y2];
        end
        BA_params.hyph_pos=hyph_pos;
        boat_pos(:,1)=GPSandPosition_table.Boat_pos_x_m;
        boat_pos(:,2)=GPSandPosition_table.Boat_pos_y_m;
        BA_params.boat_pos=boat_pos;
        BA_params.boat_start_latlong=[GPSandPosition_table.Boat_Latitude(1),GPSandPosition_table.Boat_Longitude(1)];

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////


%/////////// PARAMETERS for AMBIGUITY SURFACE computation ////////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%---------------- Ambiguity surface STD ----------------
AS_params.sig =0.003; %Determine sigma (standard deviation) for the Gaussian
%sig=0.0024; % From Yvonnes paper
AS_params.sig_hyperbolas = 0.00003; % STD for plotting intersecting hyperbolas
%(for visual assesment of how they are crossing)- needs to be small

%----------------  Parameters for modeled TDOA ---------------- 
% - Grid range and resolution:
minxrange= min(boat_pos(:,1))- 400; %- 400m to include the array behind the boat
maxxrange= max(boat_pos(:,1)) + 2000; %add 2km in front of the array
AS_params.xrange=[minxrange,maxxrange]; % x range in m [min,max]
AS_params.yrange=[-300,10000]; % y range in m (since boat track will be rotated
% around x-axis, and localization is bi-ambigous specify positive y-range 
% or in case sources are near trackline can add a few m on the other side- e.g. -300m).
% The max should be maximum distance you expect to still localize your signal 
AS_params.dx=10; % grid step size in m (resolution) in x direction
AS_params.dy=10;% grid step size in m (resolution) in y direction
% - Speed of sound
AS_params.c=parameters.c; % We are at the moment using the same as for TDOA tracking
%parameters.c.

%---------------- Whale swim speed ----------------
AS_params.mean_horiz_swimspeed= 0.5; %mean horizontal swim speed [m/s] 

%------------ Specify cutoff criteria for tracks around the beam ---------
% Specify what is the first degrees to be considered around the beam - this
% will affect which tracks you consider for localization
AS_params.bearing_cutoff = 60; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------- Compute TDOA cutoff based on the specified bearing cutoff--------
bear2tdoa= @(x) cosd(x).*(parameters.d/AS_params.c);
AS_params.tdoa_cutoff = bear2tdoa(AS_params.bearing_cutoff);

%----------------Create a grid for surface evaluation----------------
% At the moment grid is fixed, but in future it should move with sensors.
x=AS_params.xrange(1):AS_params.dx:AS_params.xrange(2);
AS_params.Ngp_x=length(x);
y=AS_params.yrange(1):AS_params.dy:AS_params.yrange(2);
AS_params.Ngp_y=length(y);
[AS_params.X,AS_params.Y] = meshgrid(x,y);
AS_params.wpos=[AS_params.X(:),AS_params.Y(:)]; %gives grid of N x [x,y,z] coordinates
AS_params.N= size(AS_params.wpos,1); %number of grid points to evaluate
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////








end