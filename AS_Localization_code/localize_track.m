function [SelectedTracks,AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table] = localize_track(Tracks,AS_params,BA_params,boat_pos,boat_start_latlong,hyph_pos,timevec)
%function localize_tracks.m localizes a selected TDOA track (or a gourp of
%track segments) using the ambiguity surfaces (AS) and returns computed surfaces
%along with the localization and perpendicular distance information.
%
%INPUTS:
% - Tracks: a structure containing TDOA tracks. Structure has 3 fields:
%                   ~ time - a vector of times (starting at 0) for a given
%                           track
%                   ~ time_local- a vector of times (serial date format) for a
%                           given track                          
%                   ~ tdoa - a vector of tdoas for a given track
% - AS_params: a structure containing parameters for ambiguity surface 
%              computation. (see specify_parameters.m for info on fields in
%              this structure)  
% - BA_params:a structure containing parameters for the boat and array.
%           (see specify_parameters.m for info on fields in this structure)  
% - boat_pos : boat position (in relative coordinates)- M x 2 matrix, 
%           where M is number of time steps, and each row is [x,y] coordinate 
%           for that time step
% - boat_start_latlong: boat start position in latitude and longitude 
%              (decimal degrees), 1 x 2 vector [latitude, longitude];
% - hyph_pos: two sensor positions per time step- a 2 x N x M array, where
%           N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps.         
% - timevec : a 1 x M vector of times that covers the duration of the encounter  
%
%OUTPUTS:
% - AStotal: final ambiguity surface- A x B matrix 
%           where A is number of y coordinates, and B is number of x coordinates. 
% - ASdilatetotal: final ambiguity surface with dilation, same dimensions 
%                   as AStotal. 
% - AStotal_hyperbolas: final surface with intersecting hyperbolas, same
%                       dimensions as AStotal.
% - Loc_table: table with localization information with 7 columns:
%              ~ 'TrackID'- selected track indices corresponding to entries
%                           of input structure 'Tracks'.
%              ~ 'Loc_m'- localization estimate in Cartesian coordinates 
%                         with respect to start of boat track (in m) 
%              ~ 'Loc_LatLong' - localization estimate as latitude &
%                               longitude (in decimal degrees)
%              ~ 'Loc_m_dilated' - localization estimate from dilated AS   
%                               in Cartesian coordinates with respect to  
%                               start of boat track (in m)
%              ~ 'Loc_LatLong_dilated'- localization estimate from dilated 
%                                       AS as latitude & longitude (in 
%                                       decimal degrees)
%              ~ 'distance_m' - perpendicular distance from localization 
%                              estimate ('Loc_m') to boat trackline
%              ~ 'distance_m_dilated' - perpendicular distance from 
%                                   localization estimate from dilated AS 
%                                   ('Loc_m_dilated') to boat trackline
%
%
%Pina Gruden January 2023




%//////////////////////////////////////////////////////////////////////////
%% /////////////////////1) Select TDOA track ///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for

[tdoa_measured_select,selected_indx,SelectedTracks] = select_tracks(Tracks, ...
    timevec,AS_params.tdoa_cutoff);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Compute Ambiguity Surface ///////////////////
tic
[AStotal,ASdilatetotal,AStotal_hyperbolas]= computeAS(tdoa_measured_select, ...
    selected_indx,hyph_pos,AS_params,BA_params);
toc
%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 3) Extract Localization info ///////////////////

% Select two largest peaks:
maxval=maxk(AStotal(:),2); %find two largest peaks in the ambiguity surface
ind_max_peak = [find(AStotal==maxval(1));find(AStotal==maxval(2))]; %indices of 2 maximum peaks
estim_location_m = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong = M2LatLon(estim_location_m,[boat_start_latlong(1),boat_start_latlong(2)]);

maxval=maxk(ASdilatetotal(:),2); %find two largest peaks in the dilated ambiguity surface
ind_max_peak = [find(ASdilatetotal==maxval(1));find(ASdilatetotal==maxval(2))]; %indices of 2 maximum peaks
estim_location_m_dilated = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong_dilated = M2LatLon(estim_location_m_dilated,[boat_start_latlong(1),boat_start_latlong(2)]);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 4) Compute perpendicular distance ///////////////////

%g et coefficients for a line equation for a boat:
p=polyfit(boat_pos(:,1),boat_pos(:,2),1);
% p =[p1,p0] - where p1= slope and p0= intercept (coeffs of p are in
% descending powers)

% compute perpendicular distance
d1=abs(p(1)*estim_location_m(:,1) - estim_location_m(:,2) + p(2))./sqrt(p(1)^2+1^2);
d2=abs(p(1)*estim_location_m_dilated(:,1) - estim_location_m_dilated(:,2) + p(2))./sqrt(p(1)^2+1^2);

% get perpendicular line equation (slope and intercept):
slope_perp =  -1/p(1); %slope of perpendicular line is equal to negative inverse of the slope of the line that it is perpendicular to.
interc_perp = estim_location_m(:,2) - slope_perp.*estim_location_m(:,1); 


%Create a table with localization info:
Loc_table = table(repmat({SelectedTracks},2,1),estim_location_m,estim_location_latlong, ...
    estim_location_m_dilated,estim_location_latlong_dilated,d1,d2, ...
    'VariableNames',{'TrackID','Loc_m','Loc_LatLong','Loc_m_dilated', ...
    'Loc_LatLong_dilated','distance_m','distance_m_dilated'});

end