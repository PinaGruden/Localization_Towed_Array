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

%% //////////////////// SET parameters //////////////////// 

%-------- Parameters for Ambiguity surface computation -------------
sig =0.003; %Determine sigma (standard deviation) for the Gaussian
%sig=0.0024; % From Yvonnes paper
sig_hyperbolas = 0.00003; % STD for plotting intersecting hyperbolas
%(for visual assesment of how they are crossing)- needs to be small
mean_horiz_swimspeed= 0.5; %mean horizontal swim speed (for sperm whale) [m/s] 

%--------  Parameters for the boat and array  ---------------- 
% Choose either: 
% hyph_pos = 1 to simulate boat movement and hydrophone positions
% hyph_pos = 2 to use the real data

get_hyph_pos = 2; 
switch get_hyph_pos
    case 1
        %d=parameters.d; %distance between sensors
        boatspeed_kts = 10; % boat speed in knots
        boatspeed= boatspeed_kts/1.944; %boat speed in m/s
        timestep = parameters.dt; % how much time elapses in each time step in s
        %specify line gradient and intersect along which the boat moves:
        mb=-0.5; % line gradient
        b=10; %intersect with y axis
        
    case 2
       GPSandPosition_table=readtable('GPSandPosition_table.csv'); 
       timestep = parameters.dt; % how much time elapses in each time step in s
end
%--------  Parameters for modeled TDOA ---------------- 
% - Grid range and resolution:
xrange=[-1000,4000]; % x range in m
yrange=[-4000,4000]; % y range in m
dx=10; % grid step size in m (resolution) in x direction
dy=dx;% grid step size in m (resolution) in y direction
% - Speed of sound
%c=1500;


%////////////////////Create a grid for surface evaluation//////////////////
% At the moment grid is fixed, but in future it should move with sensors.
x=xrange(1):dx:xrange(2);
Ngp_x=length(x);
y=yrange(1):dy:yrange(2);
Ngp_y=length(y);
[X,Y] = meshgrid(x,y);
wpos=[X(:),Y(:)]; %gives grid of N x [x,y,z] coordinates
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
% Pamguard tdoas multiplied by -1 to match the reference we use in the TDOA package
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


%% ////////// Get Hydrophone positions ///////////

switch get_hyph_pos
    case 1
        % Simulate Boat movement and Hydrophone positions
        % get boat positions along the specified line:
        x_line=0:1:Ntsteps+1;
        y_line=mb.*x_line+b;
        dxy_line=sqrt((x_line(2)-x_line(1))^2+(y_line(2)-y_line(1))^2);

        boatmove_m = boatspeed*timestep; %distance the boat traveles in one time step
        x_boat= x_line(1:end-1) + boatmove_m.*(x_line(2:end)-x_line(1:end-1))./dxy_line;
        y_boat=mb.*x_boat+b;

        hyph_pos=nan(2,2,Ntsteps);%[x1,y1; x2,y2];
        %first sensor is at the position of the boat
        hyph_pos(1,1,:)=x_boat(1:end-1);
        hyph_pos(1,2,:)=y_boat(1:end-1);
        %second sensor is distance d behind first
        %hyph_pos(2,1,:)=hyph_pos(1,1,1:end-1) - d.*(hyph_pos(1,1,2:end)-hyph_pos(1,1,1:end-1))./dxy_line;
        hyph_pos(2,1,:)=x_boat(1:end-1) - parameters.d.*(x_boat(2:end)-x_boat(1:end-1))./dxy_line;
        hyph_pos(2,2,:)=mb.*hyph_pos(2,1,:)+b;

    case 2
        % Get Hydrophone positions from data
         hyph_pos=nan(2,2,Ntsteps);
        for t=1:Ntsteps
            hyph_pos(:,:,t)= [GPSandPosition_table.Sensor1_pos_x_m(t),...
                GPSandPosition_table.Sensor1_pos_y_m(t);...
                GPSandPosition_table.Sensor2_pos_x_m(t),...
                GPSandPosition_table.Sensor2_pos_y_m(t)]; %[x1,y1; x2,y2];
        end

end


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

%% ------------- Compute the Surface --------------
%~~~~~~~~~ Get time step when animal is closest to 90deg (tdoa=0s)~~~~~~~~~
% This is for the purpose of Surface Dilation
[~,ind]=min(abs(tdoa_measured_select));
tsteps=1:1:Ntsteps;
tstep90=tsteps(ind); 

% Pre-allocate
LS_select= nan(1,N,Ntsteps);
LS_Hyperbolas = LS_select;
LSdilate = nan(Ngp_y,Ngp_x,Ntsteps); %swap x and y inpuput arguments since 
% reshape() will be used to fill in the values

tic
count=1; plotf=0;
for t=1:Ntsteps % for each time step compute LS
    if selected_indx(t)
        rp=hyph_pos(:,:,t);
        ip1=1;ip2=2;
        % for each hypothetical grid position wpos(wpi) this calculates
        % time between position and hydrophone position - dt1 (between
        % hyph 1 and position wpos(wpi)), dt2 (between hyph2 and position wpos(wpi)), then
        % tdoa_model is tdoa between the hyph1 and hyph2 if the source
        % is at wpos(wpi).
        dt1= 1/parameters.c.*sqrt((rp(ip1,1)-wpos(:,1)).^2 +(rp(ip1,2)-wpos(:,2)).^2);
        dt2= 1/parameters.c.*sqrt((rp(ip2,1)-wpos(:,1)).^2 +(rp(ip2,2)-wpos(:,2)).^2);
        tdoa_model=dt1-dt2;

        tdoa_diff=(tdoa_model-(tdoa_measured_select(:,t))).^2;
        LS_select(:,:,t)= exp(-1/(2*sig^2).*tdoa_diff);
        LS_Hyperbolas(:,:,t)= exp(-1/(2*sig_hyperbolas^2).*tdoa_diff);

        %///////////////////DILATE SURFACES////////////////////////////////
        % Apply Surface dilation filter (imdilate) to compensate for whale movement:
        LStotal_temp=reshape(LS_select(:,:,t),[Ngp_y,Ngp_x]);
        dt90 = timestep*abs(t-tstep90); %Elapsed time from when animals are at 90deg (0 tdoa) [s]
        if dt90==0

            LSdilate(:,:,t)=LStotal_temp; %do not dilate if animal is at the beam (the tdoas should be correct)

        else
            mdwh=dt90*mean_horiz_swimspeed; %maximum distance whale could have travelled horizontally [m]

            %gridpoints for filter in x (and also y since the same assumptions re swim speed and grid space)
            xgridf=0:dx:(mdwh+dx);
            filt_x_grid = [-fliplr(xgridf(2:end)),xgridf];
            ygridf=0:dy:(mdwh+dy);
            filt_y_grid = [-fliplr(ygridf(2:end)),ygridf];
            [Fx,Fy] = meshgrid(filt_x_grid,filt_y_grid);
            De = sqrt(Fx.^2+Fy.^2)/mdwh; %normalize to max distance that the whale could have swam (feasible distance will be 1 or lower)

            LSdilate(:,:,t)=imdilate(LStotal_temp,(De<=1));
        end
        %//////////////////////////////////////////////////////////////////

        if any(count==plotf)
            figure,hold on
            LStotal_temp=reshape(LS_select(:,:,t),[Ngp_y,Ngp_x]);
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
toc
%% ~~~~~~~~~~~~~~~~PLOT final Ambiguity Surface without Dilation~~~~~~~~~~~~
LStotal_temp=prod(LS_select,3,'omitnan');
LStotal=reshape(LStotal_temp,[Ngp_y,Ngp_x]);
figure; hold on;
s=pcolor(X,Y,LStotal);
s.EdgeColor='none';
clim([0,1])
colorbar
axis equal
hph=1;
plot(squeeze(hyph_pos(hph,1,:)),squeeze(hyph_pos(hph,2,:)),'r-','Linewidth', 3)
%plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-', 'Linewidth', 3),hold on
switch get_hyph_pos
    case 1
        plot(x_boat, y_boat,'g-','Linewidth', 1.5)
    case 2
        plot(GPSandPosition_table.Boat_pos_x_m, GPSandPosition_table.Boat_pos_y_m,'g-','Linewidth', 1.5)
end
xlabel(' x (m)'),ylabel('y (m)')
title (['Total ambiguity surface for source ', num2str(SelectedTracks)])
legend('Ambiguity Surface',['Sensor ', num2str(hph),' position'],'Boat track')
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
set(gca,'FontSize',16)

%~~~~~~~~~~~~~~~~PLOT final Ambiguity Surface with Dilation~~~~~~~~~~~~~
LSdilatetotal=prod(LSdilate,3,'omitnan');
h=figure; hold on,
s=pcolor(X,Y,LSdilatetotal);
s.EdgeColor='none';
clim([0,1])
colorbar
axis equal
hph=1;
plot(squeeze(hyph_pos(hph,1,:)),squeeze(hyph_pos(hph,2,:)),'r-','Linewidth', 3)
switch get_hyph_pos
    case 1
        plot(x_boat, y_boat,'g-','Linewidth', 1.5)
    case 2
        plot(GPSandPosition_table.Boat_pos_x_m, GPSandPosition_table.Boat_pos_y_m,'g-','Linewidth', 1.5)
end
xlabel(' x (m)'),ylabel('y (m)')
title(['Total ambiguity surface WITH dilation for source ', num2str(SelectedTracks)])
legend('Ambiguity Surface',['Sensor ', num2str(hph),' position'],'Boat track')
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
fontsize(h,16,'points')

%~~~~~~~~~~~~~~~~~~~PLOT Intersecting Hyperbolas~~~~~~~~~~~~~~~~~~~
LStotal_hyperbolas_temp = sum(LS_Hyperbolas,3,'omitnan');
LStotal_hyperbolas = reshape(LStotal_hyperbolas_temp,[Ngp_y,Ngp_x]);
figure;
s=pcolor(X,Y,LStotal_hyperbolas); hold on
% set(gca,'YDir', 'normal');
s.EdgeColor='none';
clim([0,1])
colorbar
axis equal
hph=1;
plot(squeeze(hyph_pos(hph,1,:)),squeeze(hyph_pos(hph,2,:)),'r-','Linewidth', 3)
%plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-', 'Linewidth', 3),hold on
switch get_hyph_pos
    case 1
        plot(x_boat, y_boat,'g-','Linewidth', 1.5)
    case 2
        plot(GPSandPosition_table.Boat_pos_x_m, GPSandPosition_table.Boat_pos_y_m,'g-','Linewidth', 1.5)
end
xlabel(' x (m)'),ylabel('y (m)')
legend('Intesecting hyperbolas',['Sensor ', num2str(hph),' position'], 'Boat track')
title (['Intersecting hyperbolas for source ', num2str(SelectedTracks)])
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
set(gca,'FontSize',16)