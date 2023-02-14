% Localize with Ambiguity surfaces


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had first run A0_SetUp.m!!!!! and if GPS data avaialbe also
% run A0_GetGPSData.m!

% It loads the following variables to
% Matlab workspace:
% - AS_params = a structure containing parameters for ambiguity surface 
%              computation. (see specify_parameters.m for info on fields in
%              this structure).  
% - BA_params = a structure containing parameters for the boat and array.
%           (see specify_parameters.m for info on fields in this structure)  
% - boat_pos = boat position (in relative coordinates)- M x 2 matrix, 
%           where M is number of time steps, and each row is [x,y] coordinate 
%           for that time step
% - boat_start_latlong: boat start position in latitude and longitude 
%              (decimal degrees), 1 x 2 vector [latitude, longitude];
% - hyph_pos = sensor positions per time step- a 2 x N x M array, where
%           N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps.
% - timevec = a vector of times that covers the duration of the encounter
% - Tracks_selected = a structure of tracks that satisfy the cutoff
%           criteria around the beam. It contains the following fields:
%                   ~ time - a vector of times (starting at 0) for a given
%                           track
%                   ~ time_local- a vector of times (datetime format) for a
%                           given track                          
%                   ~ tdoa - a vector of tdoas for a given track
% - folder2save2 = path to folder where results are saved to.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%//////////////////////////////////////////////////////////////////////////
%% /////////////////////1) Localize TDOA tracks ///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for, compute localization estimate and perpendicular
% distance to trackline - results are saved in 'Loc_table'

[SelectedTracks,AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table] = localizetracks(Tracks_selected, ...
    AS_params,BA_params,boat_pos,boat_start_latlong,hyph_pos,timevec);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Plot FINAL SURFACES ///////////////////

figure,
% Final ambiguity surface
subplot(131)
AStotal_S=sum(cat(3,AStotal{:}),3);
plot_AS(AStotal_S,AS_params,hyph_pos,boat_pos,Loc_table.Loc_m);
title ('Total ambiguity surface for all sources')

% Final ambiguity surface with dilation
subplot(132)
ASdilatetotal_S=sum(cat(3,ASdilatetotal{:}),3);
plot_AS(ASdilatetotal_S,AS_params,hyph_pos,boat_pos,Loc_table.Loc_m_dilated);
title ('Total ambiguity surface with dilation for all sources')

% Final surface with intersecting hyperbolas
subplot(133)
AStotal_hyperbolas_S=sum(cat(3,AStotal_hyperbolas{:}),3);
plot_AS(AStotal_hyperbolas_S,AS_params,hyph_pos,boat_pos,Loc_table.Loc_m_dilated);
title ('Intersecting hyperbolas for all sources')

%Make the figure big
ss = get(0, 'Screensize'); %
set(gcf, 'Position', [ss(1) ss(4)/2 ss(3) ss(4)/2]);


%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 3) Save ///////////////////

save([folder2save2,'Results.mat'],'Loc_table','AStotal','ASdilatetotal','AStotal_hyperbolas')
