% Localize with Ambiguity surfaces


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had first run A0_SetUp.m!!!!!

% It loads the following variables to
% Matlab workspace:
% - AS_params = a structure containing parameters for ambiguity surface 
%              computation. (see specify_parameters.m for info on fields in
%              this structure).  
% - BA_params = a structure containing parameters for the boat and array.
%           (see specify_parameters.m for info on fields in this structure)  
% - GPSandPosition_table = a table containing GPS and sensor information.
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

[AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table] = localize_tracks(Tracks_selected, ...
    AS_params,BA_params,GPSandPosition_table,hyph_pos,timevec);


% fprintf(['Estimated location for source [', num2str(SelectedTracks),'] from non-dilated AS ' ...
%     'is:\n [', repmat('%g, ', 1, numel(estim_location_m)-1), '%g] m \n and ' ...
%     '[', repmat('%f, ', 1, numel(estim_location_m)-1), '%f]  degrees \n'], ...
%     estim_location_m,estim_location_latlong)
%  
% fprintf(['Estimated location for source [', num2str(SelectedTracks),'] from dilated AS ' ...
%     'is:\n [', repmat('%g, ', 1, numel(estim_location_m_dilated)-1), '%g] m \n and ' ...
%     '[', repmat('%f, ', 1, numel(estim_location_m_dilated)-1), '%f]  degrees \n'], ...
%     estim_location_m_dilated,estim_location_latlong_dilated)
% 
% fprintf(['Estimated perpendicular distance between trackline and source [', ...
%     num2str(SelectedTracks),'] is: \n ', ...
%     '%g m and %g m for non-dilated and dilated AS estimates, respectively. \n'], ...
%     d1,d2)

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Plot ///////////////////

switch BA_params.get_hyph_pos
    case 1 % SIMUALTED GPS DATA
        boat_pos(:,1)=x_boat;
        boat_pos(:,2)=y_boat;

    case 2 % REAL GPS DATA
        boat_pos(:,1)=GPSandPosition_table.Boat_pos_x_m;
        boat_pos(:,2)=GPSandPosition_table.Boat_pos_y_m;
end

% Final ambiguity surface
for k=1:size(Loc_table,1)
plot_AS(AStotal{k},AS_params,hyph_pos,boat_pos,Loc_table.Loc_m{k});
title (['Total ambiguity surface for source ', num2str(Loc_table.TrackID{k})])

% Final ambiguity surface with dilation
plot_AS(ASdilatetotal{k},AS_params,hyph_pos,boat_pos,Loc_table.Loc_m_dilated{k});
title(['Total ambiguity surface WITH dilation for source ', num2str(Loc_table.TrackID{k})])

% Final surface with intersecting hyperbolas
plot_AS(AStotal_hyperbolas{k},AS_params,hyph_pos,boat_pos,Loc_table.Loc_m{k});
title (['Intersecting hyperbolas for source ', num2str(Loc_table.TrackID{k})])
end

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 3) Save ///////////////////

save([folder2save2,'Results.mat'],'Loc_table','AStotal','ASdilatetotal','AStotal_hyperbolas')
