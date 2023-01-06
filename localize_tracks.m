function [AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table] = localize_tracks(Tracks,AS_params,BA_params,GPSandPosition_table,hyph_pos,timevec)
%function localize_tracks.m localizes a selected TDOA track (or a gourp of
%track segments) using the ambiguity surfaces (AS) and returns computed surfaces
%along with the localization and perpendicular distance information.

%INPUTS:
% - Tracks: a structure containing TDOA tracks. Structure has 3 fields:
%                   ~ time - a vector of times (starting at 0) for a given
%                           track
%                   ~ time_local- a vector of times (datetime format) for a
%                           given track                          
%                   ~ tdoa - a vector of tdoas for a given track
% - AS_params: a structure containing parameters for ambiguity surface 
%              computation. (see specify_parameters.m for info on fields in
%              this structure)  
% - BA_params:a structure containing parameters for the boat and array.
%           (see specify_parameters.m for info on fields in this structure)  
% - GPSandPosition_table,
% - hyph_pos: sensor positions per time step- a 2 x N x M array, where
%           N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps.         
% - timevec : a 1 x M vector of times that covers the duration of the encounter  


%OUTPUTS:
% - AStotal: final ambiguity surface, a 1 x M cell array where M is a 
%           number tracked sources/groups. Each cell is an A x B matrix 
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
%                              estimate ('Loc_m') to boat trackline (in m)
%              ~ 'distance_m_dilated' - perpendicular distance from 
%                                   localization estimate from dilated AS 
%                                   ('Loc_m_dilated') to boat trackline (in m)


%Pina Gruden January 2023



N=size(Tracks,2);
SelectedTracks=cell(N,1);
AStotal=cell(N,1);
ASdilatetotal=cell(N,1);
AStotal_hyperbolas=cell(N,1);
estim_location_m=cell(N,1);
estim_location_latlong=cell(N,1);
estim_location_m_dilated=cell(N,1);
estim_location_latlong_dilated=cell(N,1);
d1=cell(N,1);
d2=cell(N,1);

continue_localize = 1;
k=1;

while continue_localize ~=0 %to allow multipel groups to be localized consecutevly

%//////////////////////////////////////////////////////////////////////////
%% /////////////////////1) Select TDOA track ///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for

[tdoa_measured_select,selected_indx,SelectedTracks{k}] = select_tracks(Tracks, ...
    timevec,AS_params.tdoa_cutoff);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Compute Ambiguity Surface ///////////////////

[AStotal{k},ASdilatetotal{k},AStotal_hyperbolas{k}]= computeAS(tdoa_measured_select, ...
    selected_indx,hyph_pos,AS_params,BA_params);


%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 3) Extract Localization info ///////////////////

% NEED TO CHANGE THIS SECTION- at the moment I'm just selecting one max
% peak... In theory the peak on the other side should be at the same
% location but it's not in reality..

% % ONE WAY:
% [peakdata] =findpeaks2D(AS_params.X,AS_params.Y,AStotal);
% 
% disp(' ')
% disp(['Estimated location for source [', num2str(SelectedTracks),'] is:'])
% disp(' ')
% 
% %For the moment take just one of the estimated peaks- CHANGE THIS
% [~,n]=max(peakdata.peakZ);
% estimated_position=[peakdata.peakX(n),peakdata.peakY(n)]; %were at the moment considering 2D, so z coordinate = 0.
% fprintf(['Estimated source location is: [', repmat('%g, ', 1, numel(estimated_position)-1), '%g]\n'],estimated_position)
% fprintf('Width of the peak in X dirextion %.2f \n',peakdata.peakXWidth(n))
% fprintf('Width of the peak in Y dirextion %.2f \n',peakdata.peakYWidth(n))
% fprintf('Ambiguity value for the estimated location is %.2f \n',peakdata.peakZ(n))
% hold on, plot3(peakdata.peakX(n),peakdata.peakY(n),peakdata.peakZ(n), 'k*', 'MarkerSize',6)

% ANOTHER WAY:
ind_max_peak = find(AStotal{k} == max(max(AStotal{k}))); %corresponds to max of x y z. turns into column
estim_location_m{k} = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong{k} = M2LatLon(estim_location_m{k},[GPSandPosition_table.Boat_Latitude(1),GPSandPosition_table.Boat_Longitude(1)]);

ind_max_peak = find(ASdilatetotal{k} == max(max(ASdilatetotal{k}))); %corresponds to max of x y z. turns into column
estim_location_m_dilated{k} = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong_dilated{k} = M2LatLon(estim_location_m_dilated{k},[GPSandPosition_table.Boat_Latitude(1),GPSandPosition_table.Boat_Longitude(1)]);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 4) Compute perpendicular distance ///////////////////

%get coefficients for a line equation for a boat:
p=polyfit(GPSandPosition_table.Boat_pos_x_m,GPSandPosition_table.Boat_pos_y_m,1);
% p =[p1,p0] - where p1= slope and p0= intercept (coeffs of p are in
% descending powers)

d1{k}=abs(p(1)*estim_location_m{k}(1) - estim_location_m{k}(2) + p(2))./sqrt(p(1)^2+1^2);
d2{k}=abs(p(1)*estim_location_m_dilated{k}(1) - estim_location_m_dilated{k}(2) + p(2))./sqrt(p(1)^2+1^2);

%//////////////////////////////////////////////////////////////////////////

    continue_localize = input(['Localize another group? ' ...
        '(1 for yes, 0 for no), then enter ']);
    if continue_localize~=0 && continue_localize~=1
        error(['Not a valid choice. Select 0 to stop localizing sources,' ...
            ' or select 1 to localize another group'])
    end
    if continue_localize==1
    k=k+1;
    end

end
%Prune away all unnecessary cells:
SelectedTracks=SelectedTracks(1:k);
AStotal=AStotal(1:k);
ASdilatetotal=ASdilatetotal(1:k);
AStotal_hyperbolas=AStotal_hyperbolas(1:k);
estim_location_m=estim_location_m(1:k);
estim_location_latlong=estim_location_latlong(1:k);
estim_location_m_dilated=estim_location_m_dilated(1:k);
estim_location_latlong_dilated=estim_location_latlong_dilated(1:k);
d1=d1(1:k);
d2=d2(1:k);

%Create a table with localization info:
Loc_table = table(SelectedTracks,estim_location_m,estim_location_latlong, ...
    estim_location_m_dilated,estim_location_latlong_dilated,d1,d2, ...
    'VariableNames',{'TrackID','Loc_m','Loc_LatLong','Loc_m_dilated', ...
    'Loc_LatLong_dilated','distance_m','distance_m_dilated'});

end