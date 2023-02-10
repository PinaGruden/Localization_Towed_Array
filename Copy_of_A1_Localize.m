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

% Get boat positions
switch BA_params.get_hyph_pos
    case 1 % SIMUALTED GPS DATA
        boat_pos(:,1)=x_boat;
        boat_pos(:,2)=y_boat;

    case 2 % REAL GPS DATA
        boat_pos(:,1)=GPSandPosition_table.Boat_pos_x_m;
        boat_pos(:,2)=GPSandPosition_table.Boat_pos_y_m;
end

% Pre-allocate
N=size(Tracks_selected,2);
SelectedTracks=cell(N,1);
AStotal=cell(N,1);
ASdilatetotal=cell(N,1);
AStotal_hyperbolas=cell(N,1);
Loc_table=[];
r_col=[0.8500 0.3250 0.0980];
g_col=[0,0,0]+0.6;
y_col=[0.9290 0.6940 0.1250];

continue_localize = 1;
k=1;

while continue_localize ~=0 %to allow multiple groups to be localized consecutevly

[SelectedTracks{k},AStotal{k},ASdilatetotal{k},AStotal_hyperbolas{k},Loc_table_temp] = localize_tracks1(Tracks_selected, ...
    AS_params,BA_params,GPSandPosition_table,hyph_pos,timevec);

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
plot_AS(AStotal_hyperbolas{k},AS_params,hyph_pos,boat_pos,Loc_table_temp.Loc_m);
title (['Intersecting hyperbolas for source ', num2str(Loc_table_temp.TrackID{:})])

subplot(132)
% Final ambiguity surface
plot_AS(AStotal{k},AS_params,hyph_pos,boat_pos,Loc_table_temp.Loc_m);
title (['Total ambiguity surface for source ', num2str(Loc_table_temp.TrackID{:})])

subplot(133)
% Final ambiguity surface with dilation
plot_AS(ASdilatetotal{k},AS_params,hyph_pos,boat_pos,Loc_table_temp.Loc_m_dilated);
title(['Total ambiguity surface WITH dilation for source ', num2str(Loc_table_temp.TrackID{:})])

%Make the figure big
ss = get(0, 'Screensize'); %
set(gcf, 'Position', [ss(1) ss(4)/2 ss(3) ss(4)/2]);

drawnow
% /////////////// Keep the selected group and localizations? ///////////

keep_loc = input(['Do you want to keep this localization? ' ...
        '(1 for yes, 0 for no), then enter ']);
    if keep_loc~=0 && keep_loc~=1
        error(['Not a valid choice. Select 1 to keep localizization, or' ...
            ' or select 0 to discard'])
    end
    if keep_loc==1
        Loc_table = [Loc_table;Loc_table_temp];
    else
        % REMOVE values {k} from AS
        SelectedTracks{k}=[];
        AStotal{k}=[];
        ASdilatetotal{k}=[];
        AStotal_hyperbolas{k}=[];

        %close the plotted localizations
        close(fig) 

        % plot over selected tdoa tracks
        figure(1), hold on
        for n=1:size(SelectedTracks{k},2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local,Tracks_selected(SelectedTracks{k}(n)).tdoa,'_','Color',y_col,'LineWidth',2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local(1),Tracks_selected(SelectedTracks{k}(n)).tdoa(1),'o','Color','k','LineWidth',1.5,'MarkerFaceColor',y)
        end

        %reset group counter
        k=k-1;
    end


%///////////////// Continue to localize? ///////////////

    continue_localize = input(['Localize another group? ' ...
        '(1 for yes, 0 for no), then enter ']);
    if continue_localize~=0 && continue_localize~=1
        error(['Not a valid choice. Select 0 to stop localizing sources,' ...
            ' or select 1 to localize another group'])
    end
    if continue_localize==1
        k=k+1;

        % make the already localized tdoa tracks a different color:
        figure(1), hold on
        for n=1:size(SelectedTracks{k},2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local,Tracks_selected(SelectedTracks{k}(n)).tdoa,'_','Color',g_col,'LineWidth',2)
            plot(Tracks_selected(SelectedTracks{k}(n)).time_local(1),Tracks_selected(SelectedTracks{k}(n)).tdoa(1),'o','Color',g_col,'LineWidth',1.5,'MarkerFaceColor',g_col)
        end
    end

end

%Prune away all unnecessary cells:
SelectedTracks=SelectedTracks(1:k);
AStotal=AStotal(1:k);
ASdilatetotal=ASdilatetotal(1:k);
AStotal_hyperbolas=AStotal_hyperbolas(1:k);

%//////////////////////////////////////////////////////////////////////////
%% //////////////////// 2) Plot FINAL SURFACES ///////////////////

figure,
% Final ambiguity surface
subplot(131)
AStotal_S=sum(cat(3,AStotal{:}),3);
plot_AS(AStotal_S,AS_params,hyph_pos,boat_pos,Loc_table.Loc_m_dilated);
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
