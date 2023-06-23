function [SelectedTracks,AStotal,ASdilatetotal,AStotal_hyperbolas,Loc_table,NewGrid] = localizetracks(Tracks_selected,AS_params, ...
    BA_params,timevec)
% localizetracks.m is a function that localizes user selected TDOA
% tracks/track segments. It plots them as they are localized and user has
% the option to either keep localizations or not. User also has an option
% to localize mulitple tracks.
%
% INPUTS:
% - Tracks_selected = a structure of tracks that satisfy the cutoff
%           criteria around the beam. It contains the following fields:
%         ~ time - a vector of times (starting at 0) for a given track                          
%         ~ time_local- a vector of times (serial date format) for a given track  
%         ~ tdoa - a vector of tdoas for a given track
% - AS_params = a structure containing parameters for ambiguity surface 
%              computation. (see specify_parameters.m for info on fields in
%              this structure).  
% - BA_params:a structure containing parameters for the boat and array 
%           (see specify_parameters.m for detailed info on fields in this 
%           structure) - it includes:
%       ~ boat_pos : boat position (in relative coordinates)- M x 2 matrix, 
%         where M is number of time steps, and each row is [x,y] coordinate 
%         for that time step
%       ~ boat_start_latlong: boat start position in latitude and longitude 
%              (decimal degrees), 1 x 2 vector [latitude, longitude];
%       ~ hyph_pos: two sensor positions per time step- a 2 x N x M array, 
%           where N=2 if (x,y) coordinates are considered or N=3 if (x,y,z)
%            coordinates are considered. M = number of time steps. 
% - timevec = a vector of times that covers the duration of the encounter
%
% OUPUTS:
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
%                              estimate ('Loc_m') to boat trackline
%              ~ 'distance_m_dilated' - perpendicular distance from 
%                                   localization estimate from dilated AS 
%                                   ('Loc_m_dilated') to boat trackline
%
%
%Pina Gruden February 2023


boat_pos=BA_params.boat_pos;
hyph_pos=BA_params.hyph_pos;

% Pre-allocate
N=size(Tracks_selected,2);
SelectedTracks=cell(N,1);
AStotal=cell(N,1);
ASdilatetotal=cell(N,1);
AStotal_hyperbolas=cell(N,1);
NewGrid=cell(N,1);
Loc_table=[];
r_col=[0.8500 0.3250 0.0980];
g_col=[0,0,0]+0.6;
y_col=[0.9290 0.6940 0.1250];

continue_localize = 1;
k=1;

while continue_localize ~=0 %to allow multiple groups to be localized consecutevly

[SelectedTracks{k},AStotal{k},ASdilatetotal{k},AStotal_hyperbolas{k},Loc_table_temp,NewGrid{k}] = localize_track(Tracks_selected, ...
    AS_params,BA_params,timevec);

% /////////////// Plot selected TDOAs and Localizations///////////

% Plot selected TDOA tracks that are being localized
figure(1), hold on
lgd = legend;
lgd.AutoUpdate = 'off';
for n=1:size(SelectedTracks{k},2)
plot(Tracks_selected(SelectedTracks{k}(n)).time_local,Tracks_selected(SelectedTracks{k}(n)).tdoa,'-','Color',r_col,'LineWidth',2)
plot(Tracks_selected(SelectedTracks{k}(n)).time_local(1),Tracks_selected(SelectedTracks{k}(n)).tdoa(1),'ro','Color',r_col,'LineWidth',1.5,'MarkerFaceColor','r')
end

% Plot ambiguity surfaces:
fig=figure;
subplot(131)
% Final surface with intersecting hyperbolas
plot_AS(AStotal_hyperbolas{k},NewGrid{k},hyph_pos,boat_pos,Loc_table_temp.Loc_m);
title (['Intersecting hyperbolas for source ', num2str(Loc_table_temp.TrackID{:})])

subplot(132)
% Final ambiguity surface
est_loc_bounds=cell(1,2);
est_loc_bounds{1}=Loc_table_temp.Loc_m_min;est_loc_bounds{2}=Loc_table_temp.Loc_m_max;
plot_AS(AStotal{k},NewGrid{k},hyph_pos,boat_pos,Loc_table_temp.Loc_m,est_loc_bounds);
title (['Total ambiguity surface for source ', num2str(Loc_table_temp.TrackID{:})])

subplot(133)
% Final ambiguity surface with dilation
est_loc_bounds=cell(1,2);
est_loc_bounds{1}=Loc_table_temp.Loc_m_dilated_min;est_loc_bounds{2}=Loc_table_temp.Loc_m_dilated_max;
plot_AS(ASdilatetotal{k},NewGrid{k},hyph_pos,boat_pos,Loc_table_temp.Loc_m_dilated,est_loc_bounds);  
title(['Total ambiguity surface WITH dilation for source ', num2str(Loc_table_temp.TrackID{:})])

%Make the figure big
ss = get(0, 'Screensize'); %
set(gcf, 'Position', [ss(1) ss(4)/2-100 ss(3) ss(4)/2]); % [left bottom width height] - move the bottom down a bit so that it's not all the way to the top.

drawnow
% /////////////// Keep the selected group and localizations? ///////////

keep_loc = input(['Do you want to keep this localization? ' ...
        '(1 for yes, 0 for no), then enter ']);
    if keep_loc~=0 && keep_loc~=1
        errordlg(['Not a valid choice. Select 1 to keep localizization, or' ...
            ' or select 0 to discard'],'Error')
        keep_loc = input(['Do you want to keep this localization? ' ...
        '(1 for yes, 0 for no), then enter ']);
    end
    if keep_loc==1
        Loc_table = [Loc_table;Loc_table_temp];
    else
        

        % plot over selected tdoa tracks
        figure(1), hold on
        for n=1:size(SelectedTracks{k},2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local,Tracks_selected(SelectedTracks{k}(n)).tdoa,'-','Color',y_col,'LineWidth',2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local(1),Tracks_selected(SelectedTracks{k}(n)).tdoa(1),'o','Color','k','LineWidth',1.5,'MarkerFaceColor','y')
        end
        drawnow

        % REMOVE values {k} from AS
        SelectedTracks{k}=[];
        AStotal{k}=[];
        ASdilatetotal{k}=[];
        AStotal_hyperbolas{k}=[];
        NewGrid{k}=[];

        %close the plotted localizations
        close(fig) 

        %reset group counter
        k=k-1;
    end


%///////////////// Continue to localize? ///////////////

    continue_localize = input(['Localize another group? ' ...
        '(1 for yes, 0 for no), then enter ']);
    if continue_localize~=0 && continue_localize~=1
        errordlg(['Not a valid choice. Select 0 to stop localizing sources,' ...
            ' or select 1 to localize another group'], 'Error')
        continue_localize = input(['Localize another group? ' ...
        '(1 for yes, 0 for no), then enter ']);
    end
    if continue_localize==1
        % make the already localized tdoa tracks a different color:
        figure(1), hold on
        for n=1:size(SelectedTracks{k},2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local,Tracks_selected(SelectedTracks{k}(n)).tdoa,'-','Color',g_col,'LineWidth',2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local(1),Tracks_selected(SelectedTracks{k}(n)).tdoa(1),'o','Color',g_col,'LineWidth',1.5,'MarkerFaceColor',g_col)
        end
        drawnow

        %increase group counter
        k=k+1;
    end

end

%Prune away all unnecessary cells:
SelectedTracks=SelectedTracks(1:k);
AStotal=AStotal(1:k);
ASdilatetotal=ASdilatetotal(1:k);
AStotal_hyperbolas=AStotal_hyperbolas(1:k);
NewGrid=NewGrid(1:k);
end