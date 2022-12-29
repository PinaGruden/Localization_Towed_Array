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


%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 4) Plot ///////////////////

switch BA_params.get_hyph_pos
    case 1 % SIMUALTED GPS DATA
        boat_pos(:,1)=x_boat;
        boat_pos(:,2)=y_boat;

    case 2 % REAL GPS DATA
        boat_pos(:,1)=GPSandPosition_table.Boat_pos_x_m;
        boat_pos(:,2)=GPSandPosition_table.Boat_pos_y_m;
end

% Final ambiguity surface
plot_AS(AStotal,AS_params,hyph_pos,boat_pos);
title (['Total ambiguity surface for source ', num2str(SelectedTracks)])

% Final ambiguity surface with dilation
plot_AS(ASdilatetotal,AS_params,hyph_pos,boat_pos);
title(['Total ambiguity surface WITH dilation for source ', num2str(SelectedTracks)])

% Final surface with intersecting hyperbolas
plot_AS(AStotal_hyperbolas,AS_params,hyph_pos,boat_pos);
title (['Intersecting hyperbolas for source ', num2str(SelectedTracks)])
