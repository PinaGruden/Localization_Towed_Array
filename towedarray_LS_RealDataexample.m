% Towed Array Source Localization with Ambiguity Surface method

% REAL DATA example:
% This example uses TDOA tracks extracted by TDOA tracking package and
% selects tracks that are within certain degrees of the beam. For those
% tracks ambiguity surfaces are computed and localizations made. Note, boat
% track and hydrophone positions are simulated.

% Make sure you have first run TDOA tracking package ("TDOA_tracking_master") to extract TDOA tracks
% Also make sure you have extracted data from Pamguard sql database with
% "Extract_Pamguard_detections" package

clear,close all

% /////////////////// ADD PATHS and EXTRACTED TRACKS: ////////////////////
% Add path for plotting- change to reflect your path to the TDOA tracking package:
addpath('/Users/pinagruden/Dropbox/Pina/Git_projects/TDOA_tracking_master'); 

% Load all Pamguard tdoas from detected whistles and clicks
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Ground_truth_fromJenn/Lasker_AC109/Lasker_AC109_Extracted_AnnotatedClicks.mat')
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Ground_truth_fromJenn/Lasker_AC109/Lasker_AC109_Extracted_AnnotatedWhistles.mat')
pamguard_parameters=parameters;

% Load data/tracks extracted by the TDOA tracking package:
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/Tracking_Towed_Array/LocalizationTest/LaskerAC109_Results/Lasker_AC109_Results.mat')
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/Tracking_Towed_Array/LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_clicks_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL')
Rxy_envelope_both{1}=Rxy_envelope_ALL.*scalar_clicks;
load('/Users/pinagruden/Dropbox/Pina/HAWAII/MATLAB/Code/Tracking_Towed_Array/LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_whistles_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL','lags','t_serialdate')
Rxy_envelope_both{2}=Rxy_envelope_ALL.*scalar_whistles;

% Check that the TDOA tracking and Pamguard used the same sensor
% separation:
if ~isequal(parameters.d,pamguard_parameters.d)
    error(['The sensor spacing between TDOA tracking and Pamguard detections' ...
        ' does not match. Re-check sensor separation and compute TDOAs for' ...
        ' Pamguard for the sensor separation that matches parameters.d.'])
end

%//////////////////// SET parameters //////////////////// 

%-------- Parameters for Ambiguity surface computation -------------
sig =0.008; %Determine sigma (standard deviation) for the Gaussian
%sig=0.0024; % From Yvonnes paper
sig_hyperbolas = 0.0003; % STD for plotting intersecting hyperbolas
%(for visual assesment of how they are crossing)- needs to be small

%--------  Parameters for the boat and array simulation ----------------  
d=parameters.d; %distance between sensors
boatspeed_kts = 10; % boat speed in knots
boatspeed= boatspeed_kts/1.944; %boat speed in m/s
timestep = parameters.dt; % how much time elapses in each time step in s

%--------  Parameters for modeled TDOA ---------------- 
% - Grid range and resolution:
xrange=[0,4000]; % x range in m
yrange=[0,4000]; % y range in m
dx=10; % grid step size in m (resolution) in x direction
dy=dx;% grid step size in m (resolution) in y direction
% - Speed of sound
c=1500;


%////////////////////Create a grid for surface evaluation//////////////////
% At the moment grid is fixed, but in future it should move with sensors.
x=xrange(1):dx:xrange(2);
Ngp_x=length(x);
y=yrange(1):dy:yrange(2);
Ngp_y=length(y);
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate


%///////////////////Select Measurements/////////////////////////
% ////// Find all tracks that cross the beam /////////////

%---------1. Search for all tracks within x degrees of the beam------------
bear2tdoa = @(x) cosd(x).*(parameters.d/parameters.c);
maxtdoa=bear2tdoa(60); %Search for +/-30 deg of the beam

indx=nan(size(Tracks,2),1);
for k=1:size(Tracks,2)
indx(k)=ismembertol(0,Tracks(k).tdoa,maxtdoa); %Checked and works good
end
indx=logical(indx);
Nsources=sum(indx);
Tracks_selected=Tracks(indx);

%---2. Identify min, max times for these tracks & create a time vector---
mint=min(vertcat(Tracks_selected(:).time));
maxt=max(vertcat(Tracks_selected(:).time));
timevec=mint:parameters.dt:maxt;
Ntsteps=numel(timevec); %number of time steps


%---------------3. PLOT All tracks and selected:------------------
% a) Plot against cross-correlogram
parameters.saveplotsofresults=0;
measure.T=0;
plot_results(t_serialdate,lags,Rxy_envelope_both,measure,Tracks, parameters)
figure(2), hold on
for k=1:Nsources
text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1),num2str(k),'Color','r')
end

% b) Plot against all Pamguard detections
figure,hold on;
plot(datenum(All_data_w.time_UTC),-1.*All_data_w.tdoa,'o', ...
    'Color',[0,0,0]+0.85, 'MarkerFaceColor',[0,0,0]+0.85, 'MarkerSize', 3)
plot(datenum(All_data_c.time_UTC),-1.*All_data_c.tdoa,'.','Color',[0,0,0]+0.85)
for k=1:size(Tracks,2)
plot(Tracks(k).time_local, Tracks(k).tdoa,'-','LineWidth',4), hold on
end
set(gca,'YDir', 'reverse');
datetick('x','keeplimits');
xlim([t_serialdate(1),t_serialdate(end)])
for k=1:Nsources
text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1),num2str(k),'Color','r')
end

%---------- 4. Re-arrange tracks in a matrix (Nsources x Ntsteps)--------
tdoa_measured= nan(Nsources,Ntsteps);
for k=1:Nsources
    time_indx=ismember(timevec,Tracks_selected(k).time);
    tdoa_measured(k,time_indx)= Tracks_selected(k).tdoa;
end


% ////////// Simulate Boat movement and Hydrophone positions ///////////

hyph_pos(:,:,1)=[0,0,0;d,0,0]; %[x1,y1,z1; x2,y2,z2];
boatmove = boatspeed*timestep;
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[boatmove,0,0;boatmove,0,0];
end
Nsensors=size(hyph_pos,1);

%% /////////////////////Compute Ambiguity Surfaces///////////////////

%-------- Select which tracks you want to compute the surfaces for------
prompt = "Which track/tracks you want to compute the ambiguity surface for? \n" + ...
    " Note, enter a single track number or if track is fragmented, \n " + ...
    "then enter fragments as a 1xN vector. \n " + ...
    "Look at Figure 2 and Figure 4 for track numbers. \n ";

SelectedTracks =input(prompt);
selected_indx=zeros(length(SelectedTracks),length(timevec));
for n=1:length(SelectedTracks)
ind=find(abs(Tracks_selected(SelectedTracks(n)).tdoa)<maxtdoa);
selected_indx(n,:)=ismember(timevec,Tracks_selected(SelectedTracks(n)).time(ind));
end
selected_indx=sum(selected_indx,1);
tdoa_measured_select=sum(tdoa_measured(SelectedTracks,:),1,'omitnan');
tdoa_measured_select(all(isnan(tdoa_measured(SelectedTracks,:)),1)) = NaN;

%------------- Compute the Surface --------------
% Pre-allocate
LS_select= nan(1,N,Ntsteps);

count=1; plotf=0;
for t=1:Ntsteps % for each time step compute LS
    if selected_indx(t)
        rp=hyph_pos(:,:,t);
        ip1=1;ip2=2;
        tdoa_model=zeros(1,N);
        for wpi=1:N % for each hypothetical grid position wpos(wpi) this calculates
            % time between position and hydrophone position - dt1 (between
            % hyph 1 and position wpos(wpi)), dt2 (between hyph2 and position wpos(wpi)), then
            % tdoa_model is tdoa between the hyph1 and hyph2 if the source
            % is at wpos(wpi).
            dt1 = 1/c*sqrt( (rp(ip1,1)-wpos(wpi,1))^2 + ...
                (rp(ip1,2)-wpos(wpi,2))^2 + ...
                (rp(ip1,3)-wpos(wpi,3))^2);
            dt2 = 1/c*sqrt( (rp(ip2,1)-wpos(wpi,1))^2 + ...
                (rp(ip2,2)-wpos(wpi,2))^2 + ...
                (rp(ip2,3)-wpos(wpi,3))^2);
            tdoa_model(wpi) = dt1-dt2;
        end

        LS_select(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-(-1.*tdoa_measured_select(:,t))).^2);

        if any(count==plotf)
            figure,hold on
            LStotal_temp=reshape(LS_select(:,:,t),[Ngp_x,Ngp_y,Ngp_z]);
            s=pcolor(X,Y,LStotal_temp);
            s.EdgeColor='none';
            clim([0,1])
            axis equal
            plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
            plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
            colorbar
            xlabel(' x (m)'),ylabel('y (m)')
            title(['LS for selected source, time step ', num2str(t)])
        end

        count=count+1;
    end
end

% ~~~~~~~~~~~~~~~~PLOT final Ambiguity Surface without Dilation~~~~~~~~~~~~
LStotal_temp=prod(LS_select,3,'omitnan');
LStotal=reshape(LStotal_temp,[Ngp_x,Ngp_x,Ngp_z]);
figure; hold on;
s=pcolor(X,Y,LStotal);
s.EdgeColor='none';
clim([0,1])
colorbar
axis equal
hph=1;
plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-', 'Linewidth', 3),hold on
xlabel(' x (m)'),ylabel('y (m)')
title (['Total ambiguity surface for source ', num2str(SelectedTracks)])
legend('Ambiguity Surface','Boat track')
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
set(gca,'FontSize',16)
