% Localize with Ambiguity surfaces


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had first run A0_SetUp.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%//////////////////////////////////////////////////////////////////////////
%% /////////////////////1) Select TDOA track ///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for

[tdoa_measured_select,selected_indx,SelectedTracks] = select_tracks(Tracks_selected, ...
    timevec,AS_params.tdoa_cutoff);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Compute Ambiguity Surface ///////////////////

[AStotal,ASdilatetotal,AStotal_hyperbolas]= computeAS(tdoa_measured_select, ...
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
ind_max_peak = find(AStotal == max(max(AStotal))); %corresponds to max of x y z. turns into column
estim_location_m = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong = M2LatLon(estim_location_m,[GPSandPosition_table.Boat_Latitude(1),GPSandPosition_table.Boat_Longitude(1)]);

fprintf(['Estimated location for source [', num2str(SelectedTracks),'] from non-dilated AS ' ...
    'is:\n [', repmat('%g, ', 1, numel(estim_location_m)-1), '%g] m \n and ' ...
    '[', repmat('%f, ', 1, numel(estim_location_m)-1), '%f]  degrees \n'], ...
    estim_location_m,estim_location_latlong)

ind_max_peak = find(ASdilatetotal == max(max(ASdilatetotal))); %corresponds to max of x y z. turns into column
estim_location_m_dilated = [AS_params.X(ind_max_peak),AS_params.Y(ind_max_peak)];
estim_location_latlong_dilated = M2LatLon(estim_location_m_dilated,[GPSandPosition_table.Boat_Latitude(1),GPSandPosition_table.Boat_Longitude(1)]);

fprintf(['Estimated location for source [', num2str(SelectedTracks),'] from dilated AS ' ...
    'is:\n [', repmat('%g, ', 1, numel(estim_location_m_dilated)-1), '%g] m \n and ' ...
    '[', repmat('%f, ', 1, numel(estim_location_m_dilated)-1), '%f]  degrees \n'], ...
    estim_location_m_dilated,estim_location_latlong_dilated)

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 4) Compute perpendicular distance ///////////////////

%get coefficients for a line equation for a boat:
p=polyfit(GPSandPosition_table.Boat_pos_x_m,GPSandPosition_table.Boat_pos_y_m,1);
% p =[p1,p0] - where p1= slope and p0= intercept (coeffs of p are in
% descending powers)

d1=abs(p(1)*estim_location_m(1) - estim_location_m(2) + p(2))./sqrt(p(1)^2+1^2);
d2=abs(p(1)*estim_location_m_dilated(1) - estim_location_m_dilated(2) + p(2))./sqrt(p(1)^2+1^2);

fprintf(['Estimated perpendicular distance between trackline and source [', ...
    num2str(SelectedTracks),'] is: \n ', ...
    '%g m and %g m for non-dilated and dilated AS estimates, respectively. \n'], ...
    d1,d2)

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 5) Plot ///////////////////

switch BA_params.get_hyph_pos
    case 1 % SIMUALTED GPS DATA
        boat_pos(:,1)=x_boat;
        boat_pos(:,2)=y_boat;

    case 2 % REAL GPS DATA
        boat_pos(:,1)=GPSandPosition_table.Boat_pos_x_m;
        boat_pos(:,2)=GPSandPosition_table.Boat_pos_y_m;
end

% Final ambiguity surface
plot_AS(AStotal,AS_params,hyph_pos,boat_pos,estim_location_m);
title (['Total ambiguity surface for source ', num2str(SelectedTracks)])

% Final ambiguity surface with dilation
plot_AS(ASdilatetotal,AS_params,hyph_pos,boat_pos,estim_location_m_dilated);
title(['Total ambiguity surface WITH dilation for source ', num2str(SelectedTracks)])

% Final surface with intersecting hyperbolas
plot_AS(AStotal_hyperbolas,AS_params,hyph_pos,boat_pos,estim_location_m);
title (['Intersecting hyperbolas for source ', num2str(SelectedTracks)])
