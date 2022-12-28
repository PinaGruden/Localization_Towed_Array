function [AS_params,BA_params] = specify_parameters(parameters)
% Specify parameters for towed array localization:
% Make sure you change/adjust parameters under captions "Changable"

%% //////////////////// SET PARAMETERS //////////////////// 

%/////////// PARAMETERS for AMBIGUITY SURFACE computation ////////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
%---------------- Ambiguity surface STD ----------------
AS_params.sig =0.003; %Determine sigma (standard deviation) for the Gaussian
%sig=0.0024; % From Yvonnes paper
AS_params.sig_hyperbolas = 0.00003; % STD for plotting intersecting hyperbolas
%(for visual assesment of how they are crossing)- needs to be small

%----------------  Parameters for modeled TDOA ---------------- 
% - Grid range and resolution:
AS_params.xrange=[-1000,4000]; % x range in m
AS_params.yrange=[-4000,4000]; % y range in m
AS_params.dx=10; % grid step size in m (resolution) in x direction
AS_params.dy=10;% grid step size in m (resolution) in y direction
% - Speed of sound
%c=1500; % We are at the moment using the same as for TDOA tracking
%parameters.c.

%---------------- Whale swim speed ----------------
AS_params.mean_horiz_swimspeed= 0.5; %mean horizontal swim speed (for sperm whale) [m/s] 

%------------ Specify cutoff criteria for tracks around the beam ---------
% Specify what is the first degrees to be considered around the beam - this
% will affect which tracks you consider for localization
AS_params.bearing_cuttof = 60;
bear2tdoa= @(x) cosd(x).*(parameters.d/parameters.c);
AS_params.tdoa_cutoff = bear2tdoa(AS_params.bearing_cuttof);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%////////////////////Create a grid for surface evaluation//////////////////
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


%///////////////////// PARAMETERS for the BOAT and ARRAY /////////////////

% ~~~~~~~~~~~~~~~~~~~~~~~CHANGABLE:~~~~~~~~~~~~~~~~~~~~~~~~
% Choose either: 
% hyph_pos = 1 to simulate boat movement and hydrophone positions.
% hyph_pos = 2 to use the real data

BA_params.get_hyph_pos = 2; 

if BA_params.get_hyph_pos==1
    boatspeed_kts = 10; % boat speed in knots
    %specify line gradient and intersect along which the boat moves:
    BA_params.mb=-0.5; % line gradient
    BA_params.b=10; %intersect with y axis
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMPUTED:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
switch BA_params.get_hyph_pos
    case 1
        BA_params.boatspeed= boatspeed_kts/1.944; %boat speed in m/s
        BA_params.timestep = parameters.dt; % how much time elapses in each time step in s 
        BA_params.d = parameters.d; % distance between the sensors
    case 2
%        BA_params.GPSandPosition_table=readtable('GPSandPosition_table.csv'); 
       BA_params.timestep = parameters.dt; % how much time elapses in each time step in s
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%///////////////////////////////////////////////////////////////////////






end