function [SelectedTracks,AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table,NewGrid] = localize_track(Tracks,AS_params,BA_params,boat_pos,boat_start_latlong,hyph_pos,timevec)
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
loc_pos=0;
while loc_pos~=1
%% /////////////////////1) Select TDOA track ///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for

[tdoa_measured_select,selected_indx,SelectedTracks] = select_tracks(Tracks, ...
    timevec,AS_params.tdoa_cutoff);

%//////////////////////////////////////////////////////////////////////////
%% /////////////////////2) ROTATE coordinate grid ///////////////////
%Rotate xy coordinates to align with x-axis 

%get coefficients for a line equation for a boat:
p=polyfit(boat_pos(:,1),boat_pos(:,2),1);
% p =[p1,p0] - where p1= slope and p0= intercept (coeffs of p are in
% descending powers)

phi=atan(p(1)); %the angle between the line and x-axis is inverse tan of the slope

rotmatfnc= @(x) [cos(x),-sin(x);sin(x),cos(x)]; %rotation matrix

% Remove y-offset
boat_pos_shift = [boat_pos(:,1),boat_pos(:,2)-p(2)];
hyph_pos_shift=hyph_pos;
hyph_pos_shift(:,2,:)= hyph_pos_shift(:,2,:)-p(2);

%get rotated boat and hydrophone positions:
%boat_pos_rotd=(rotmatfnc(-phi)*boat_pos_shift')';
hyph_pos_rotd= pagetranspose(pagemtimes(rotmatfnc(-phi),pagetranspose(hyph_pos_shift(:,:,:))));

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 3) Compute Ambiguity Surface ///////////////////
tic
[AStotal,ASdilatetotal,AStotal_hyperbolas]= computeAS(tdoa_measured_select, ...
    selected_indx,hyph_pos_rotd,AS_params,BA_params);
toc
%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 4) Extract Localization info ///////////////////

% Select tallest peak as localization (since it's biambigous)
% ----------------Non-dilated surfaces-------------------
% Marginalize along x-axis
%x_marg= max(AStotal,[],1);
x_marg= sum(AStotal,1);
[~,locs_x]=findpeaks(x_marg,AS_params.X(1,:),'SortStr', 'descend', 'NPeaks', 1); 

% Marginalize along y-axis
% y_marg=max(AStotal,[],2);%
y_marg=sum(AStotal,2);
[pks_y,locs_y]=findpeaks(y_marg,AS_params.Y(:,1),'SortStr', 'descend', 'NPeaks', 1);

estim_location_m_rotd = [locs_x(:),locs_y(:)];


%------------------ Dilated surfaces----------------------
% Marginalize along x-axis
x_marg_dilate= sum(ASdilatetotal,1);
[~,locs_x]=findpeaks(x_marg_dilate,AS_params.X(1,:),'SortStr', 'descend', 'NPeaks', 1);

% Marginalize along y-axis
y_marg_dilate=sum(ASdilatetotal,2);
[pks_y_dilate,locs_y]=findpeaks(y_marg_dilate,AS_params.Y(:,1), 'SortStr', 'descend', 'NPeaks', 1);

estim_location_m_dilated_rotd = [locs_x(:),locs_y(:)];

% Select two largest peaks:
% maxval=maxk(AStotal(:),2); %find two largest peaks in the ambiguity surface
% ind_max_peak = [find(AStotal==maxval(1));find(AStotal==maxval(2))]; %indices of 2 maximum peaks
% estim_location_m = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
% estim_location_latlong = M2LatLon(estim_location_m,[boat_start_latlong(1),boat_start_latlong(2)]);
% 
% maxval=maxk(ASdilatetotal(:),2); %find two largest peaks in the dilated ambiguity surface
% ind_max_peak = [find(ASdilatetotal==maxval(1));find(ASdilatetotal==maxval(2))]; %indices of 2 maximum peaks
% estim_location_m_dilated = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
% estim_location_latlong_dilated = M2LatLon(estim_location_m_dilated,[boat_start_latlong(1),boat_start_latlong(2)]);

% Check if localization can be obtained- i.e. there is sufficient change in
% bearings:
% 
if size(estim_location_m_rotd,2)<2 || size(estim_location_m_dilated_rotd,2)<2 % Not sufficient change in bearings - localization not possible
    loc_pos=0;
errorMsg = "Localization not possible. \n" + ...
    "Selected track(s) do not result in a sufficient bearing change for successfull localization. \n" + ...
    "Try again. \n";
fprintf(2,errorMsg)
else % Sufficient change in bearings - localization possible
    loc_pos=1;
end

end
%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 5) Compute perpendicular distance ///////////////////

% ----------------Non-dilated surfaces-------------------
d_m=estim_location_m_rotd(2); % distance is equal to y-coordinate (since boat travels along x-axis)

% compute error on the estimated distance (90% from the estimated
% distance)
ind_ymarg=find(y_marg>pks_y*0.9);
d_m_low=AS_params.Y(ind_ymarg(1),1);
d_m_high=AS_params.Y(ind_ymarg(end),1);

%------------------ Dilated surfaces----------------------
d_m_dilated=estim_location_m_dilated_rotd(2);

% compute error on the estimated distance (90% from the estimated
% distance)
ind_ymarg=find(y_marg_dilate>pks_y_dilate*0.9);
d_m_dilated_low=AS_params.Y(ind_ymarg(1),1);
d_m_dilated_high=AS_params.Y(ind_ymarg(end),1);

%-----------------------------------------------------------------
% The result can be checked by running the following (also if you dont
% rotate coordinates that's how you compute the distance):
% %g et coefficients for a line equation for a boat:
% p=polyfit(boat_pos_rotd(:,1),boat_pos_rotd(:,2),1);
% % p =[p1,p0] - where p1= slope and p0= intercept (coeffs of p are in
% % descending powers)
% 
% % compute perpendicular distance
% d_m=abs(p(1)*estim_location_m_rotd(:,1) - estim_location_m_rotd(:,2) + p(2))./sqrt(p(1)^2+1^2);
%
% % get perpendicular line equation (slope and intercept):
% slope_perp =  -1/p(1); %slope of perpendicular line is equal to negative inverse of the slope of the line that it is perpendicular to.
% interc_perp = estim_location_m(:,2) - slope_perp.*estim_location_m(:,1); 

%//////////////////////////////////////////////////////////////////////////
%% /////////////// 6) REVERSE rotation of the coordinate grid ////////////// 

%---------Reverse the rotation and recover y-offset:----------
% -------- plus get min, max positions based on min, max distance-------

%~~~~~~~~~~~~~~~~~~~~~~ Non-dialted ~~~~~~~~~~~~~~~~~~~~~ 
% - Position
estim_location_m_shift=(rotmatfnc(phi)*estim_location_m_rotd')';
estim_location_m=estim_location_m_shift;
estim_location_m(:,2)=estim_location_m_shift(:,2)+p(1);

% - Position min
estim_location_m_shift_low=(rotmatfnc(phi)*[estim_location_m_rotd(1),d_m_low]')';
estim_location_m_low=estim_location_m_shift_low;
estim_location_m_low(:,2)=estim_location_m_shift_low(:,2)+p(1);

% - Poisition max
estim_location_m_shift_hi=(rotmatfnc(phi)*[estim_location_m_rotd(1),d_m_high]')';
estim_location_m_hi=estim_location_m_shift_hi;
estim_location_m_hi(:,2)=estim_location_m_shift_hi(:,2)+p(1);


%~~~~~~~~~~~~~~~~~~~~~~ Dilated ~~~~~~~~~~~~~~~~~~~~~~
% - Position
estim_location_m_dilated_shift=(rotmatfnc(phi)*estim_location_m_dilated_rotd')';
estim_location_m_dilated=estim_location_m_dilated_shift;
estim_location_m_dilated(:,2)=estim_location_m_dilated_shift(:,2)+p(1);

% - Position min
estim_location_m_dilated_shift_low=(rotmatfnc(phi)*[estim_location_m_dilated_rotd(1),d_m_dilated_low]')';
estim_location_m_dilated_low=estim_location_m_dilated_shift_low;
estim_location_m_dilated_low(:,2)=estim_location_m_dilated_shift_low(:,2)+p(1);

% - Poisition max
estim_location_m_dilated_shift_hi=(rotmatfnc(phi)*[estim_location_m_dilated_rotd(1),d_m_dilated_high]')';
estim_location_m_dilated_hi=estim_location_m_dilated_shift_hi;
estim_location_m_dilated_hi(:,2)=estim_location_m_dilated_shift_hi(:,2)+p(1);



%-------------- Get lat long of the non-rotated localizations: ------------
estim_location_latlong = M2LatLon(estim_location_m,[boat_start_latlong(1),boat_start_latlong(2)]);
estim_location_latlong_dilated = M2LatLon(estim_location_m_dilated,[boat_start_latlong(1),boat_start_latlong(2)]);

%--------------- Get new gridpoints of non-rotated grid:------------------
new_xy= (rotmatfnc(phi)*[AS_params.wpos]')';
X_new=reshape(new_xy(:,1),size(AS_params.X));
Y_new=reshape(new_xy(:,2),size(AS_params.Y));
NewGrid.X=X_new;
NewGrid.Y=Y_new;

%//////////////////////////////////////////////////////////////////////////
%% /////////////// 7) Create a table with localization info ////////////// 

Loc_table = table({SelectedTracks}, ...
    estim_location_m,estim_location_m_low,estim_location_m_hi,estim_location_latlong, ...
    estim_location_m_dilated,estim_location_m_dilated_low,estim_location_m_dilated_hi,estim_location_latlong_dilated, ...
    d_m,[d_m_low,d_m_high],d_m_dilated,[d_m_dilated_low,d_m_dilated_high], ...
    'VariableNames',{'TrackID', ...
    'Loc_m','Loc_m_min','Loc_m_max','Loc_LatLong', ...
    'Loc_m_dilated','Loc_m_dilated_min','Loc_m_dilated_max','Loc_LatLong_dilated', ...
    'distance_m','distance_m_minmax','distance_m_dilated','distance_m_dilated_minmax'});

end