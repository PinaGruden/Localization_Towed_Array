% Compute Ambiguity surface for 1 source for towed array

% Assumes:
% - all TDOAs are from one source
% - stationary array
% - 2D
% - Isospeed profile


%//////////////////////////////////////////////////////////////////////////
%% SIMULATION 1: Stationary source, moving array, no noise on measurements
clear, close all
%Determine sigma (variance) for the Gaussian:
sig =0.001;

Ntsteps=4; %number of time steps

% Simulate Hydrophone positions
hyph_pos(:,:,1)=[0,0,0;40,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[50,0,0;50,0,0];
end

% True whale position
true_wpos= [100,50,0]; %[x,y,z];
Nsources=size(true_wpos,1);

% True TDOAs
true_tdoa=zeros(1,Ntsteps);
ip1=1;ip2=2;
c=1500;
for t=1:Ntsteps
    rp=hyph_pos(:,:,t);
    dt1 = 1/c*sqrt((rp(ip1,1)-true_wpos(1))^2 + ...
        (rp(ip1,2)-true_wpos(2))^2 + ...
        (rp(ip1,3)-true_wpos(3))^2);
    dt2 = 1/c*sqrt((rp(ip2,1)-true_wpos(1))^2 + ...
        (rp(ip2,2)-true_wpos(2))^2 + ...
        (rp(ip2,3)-true_wpos(3))^2);
    true_tdoa(t) = dt1-dt2;
end

tdoa_measured = true_tdoa; %Assume for the moment no noise in measurements

% Cell grids to evaluate:
% It could move with the hydrophone position or cover an entire area - what
% is better??
x=-200:10:200;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate

LS=nan(Nsources,N,Ntsteps);
fig1=figure; hold on
for t=1:Ntsteps % for each time step compute LS

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

LS(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-tdoa_measured(t)).^2);

LStotal_temp=reshape(LS(:,:,t),[Ngp_x,Ngp_x,Ngp_z]);
subplot(1,Ntsteps,t)
pcolor(X,Y,LStotal_temp); hold on
clim([0,1])
plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
plot(true_wpos(1),true_wpos(2),'r*')
colorbar
xlabel(' x (m)'),ylabel('y (m)')
title(['Ambiguity surface for time step ', num2str(t)])

end
fontsize(fig1,14,'points')

LStotal_temp=prod(LS,3);
LStotal=reshape(LStotal_temp,[Ngp_x,Ngp_x,Ngp_z]);
fig2=figure; hold on,
pcolor(X,Y,LStotal); 
clim([0,1])
plot([hyph_pos(ip1,1,1),hyph_pos(ip1,1,end)],[hyph_pos(ip1,2,1),hyph_pos(ip1,2,end)],'r-'),hold on
plot(true_wpos(1),true_wpos(2),'r*')
xlabel(' x (m)'),ylabel('y (m)')
title('Total ambiguity surface')
fontsize(fig2,14,'points')


%//////////////////////////////////////////////////////////////////////////
%% SIMULATION 2: Stationary source, moving array, noisy measurements 
clear, %close all
%Determine sigma (variance) for the Gaussian:
sig =0.003;

% measurement noise (needs to be equal or less than what we assume for sigma):
measure_noise=0.002;

Ntsteps=4; %number of time steps

% Simulate Hydrophone positions
d=40; %distance between sensors
hyph_pos(:,:,1)=[0,0,0;d,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[50,0,0;50,0,0];
end

% True whale position
true_wpos= [100,50,0]; %[x,y,z];
Nsources=size(true_wpos,1);

% True TDOAs
true_tdoa=zeros(1,Ntsteps);
ip1=1;ip2=2;
c=1500;
for t=1:Ntsteps
    rp=hyph_pos(:,:,t);
    dt1 = 1/c*sqrt((rp(ip1,1)-true_wpos(1))^2 + ...
        (rp(ip1,2)-true_wpos(2))^2 + ...
        (rp(ip1,3)-true_wpos(3))^2);
    dt2 = 1/c*sqrt((rp(ip2,1)-true_wpos(1))^2 + ...
        (rp(ip2,2)-true_wpos(2))^2 + ...
        (rp(ip2,3)-true_wpos(3))^2);
    true_tdoa(t) = dt1-dt2;
end

%Get measurements:
tdoa_measured = true_tdoa + randn(size(true_tdoa)).*measure_noise; 
tdoa2bearing= @(x) acosd(1500/d*x); %bear2tdoa = @(x) cosd(x).*(d/c);
fprintf('True bearings are: %.4f \n', tdoa2bearing(true_tdoa))
fprintf('Measured bearings are: %.4f \n', tdoa2bearing(tdoa_measured))


% Cell grids to evaluate:
% It could move with the hydrophone position or cover an entire area - what
% is better??
x=-200:10:200;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate

LS=nan(Nsources,N,Ntsteps);
for t=1:Ntsteps % for each time step compute LS

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

LS(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-tdoa_measured(t)).^2);

% LStotal_temp=reshape(LS(:,:,t),[Ngp_x,Ngp_x,Ngp_z]);
% figure,
% h=surf(X,Y,LStotal_temp); hold on
% colorbar
% z= get(h,'ZData');
% set(h,'ZData',z-10) 
% plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
% plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
% plot(true_wpos(1),true_wpos(2),'k*')
% xlabel(' x (m)'),ylabel('y (m)')
% shading interp
% view(0,90)

end


LStotal_temp=prod(LS,3);
LStotal=reshape(LStotal_temp,[Ngp_x,Ngp_x,Ngp_z]);
fig2=figure; hold on,
pcolor(X,Y,LStotal); 
clim([0,1])
colorbar
plot([hyph_pos(ip1,1,1),hyph_pos(ip1,1,end)],[hyph_pos(ip1,2,1),hyph_pos(ip1,2,end)],'r-'),hold on
plot(true_wpos(1),true_wpos(2),'r*','MarkerSize',6)
xlabel(' x (m)'),ylabel('y (m)')
title(['Total ambiguity surface for \sigma=',num2str(sig), ' and noise=', num2str(measure_noise)])
fontsize(fig2,14,'points')

%//////////////////////////////////////////////////////////////////////////
%% SIMULATION 3: Moving source, moving array, no noise on measurements

clear,close all
%Determine sigma (std) for the Gaussian:
sig =0.002; %0.005 (if whale moves faster)
sig_hyperbolas = 0.0003; % this is a very small std to then plot intersecting hyperbolas (visual assesment of how they are crossing)

Ntsteps=6; %number of time steps

% Simulate Hydrophone positions
hyph_pos(:,:,1)=[0,0,0;40,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[50,0,0;50,0,0];
end
Nsensors=size(hyph_pos,1);

% True whale position
true_wpos(1,:)= [100,10,0]; %[x,y,z];
for t=2:Ntsteps
true_wpos(t,:)= true_wpos(t-1,:) + [3,6,0];%[8.7,17.4,0];%%[8.7,17.4,0];%[2,5,0]; %[8.7,17.4,0] - to be equivalent
% to sperm whale swim speed given the time step used here
% [x,y,z];  + [10,16,0]; %for faster moving whale (almost 4kts)
end
Nsources=1;


% True TDOAs
true_tdoa=zeros(1,Ntsteps);
ip1=1;ip2=2;
c=1500;
for t=1:Ntsteps
    rp=hyph_pos(:,:,t);
    dt1 = 1/c*sqrt((rp(ip1,1)-true_wpos(t,1))^2 + ...
        (rp(ip1,2)-true_wpos(t,2))^2 + ...
        (rp(ip1,3)-true_wpos(t,3))^2);
    dt2 = 1/c*sqrt((rp(ip2,1)-true_wpos(t,1))^2 + ...
        (rp(ip2,2)-true_wpos(t,2))^2 + ...
        (rp(ip2,3)-true_wpos(t,3))^2);
    true_tdoa(t) = dt1-dt2;
end

tdoa_measured = true_tdoa; %Assume for the moment no noise in measurements

% Cell grids to evaluate:
% It could move with the hydrophone position or cover an entire area - what
% is better??
dx=5;
x=-300:dx:300;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate


% Get time step when animal is closest to 90deg (tdoa=0s).
[~,ind]=min(abs(tdoa_measured));
tsteps=1:1:Ntsteps;
tstep90=tsteps(ind); 

%Pre-allocate:
LS=nan(Nsources,N,Ntsteps); LS_Hyperbolas = LS;
LSdilate = nan(Ngp_x,Ngp_x,Ntsteps);
fig1=figure; hfig1=cell(1,Ntsteps);

for t=1:Ntsteps % for each time step compute LS

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

LS(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-tdoa_measured(t)).^2);
LS_Hyperbolas(:,:,t)= exp(-1/(2*sig_hyperbolas^2).*(tdoa_model-tdoa_measured(t)).^2);

LStotal_temp=reshape(LS(:,:,t),[Ngp_x,Ngp_x,Ngp_z]);

% PLOT
hfig1{t}=subplot(1,Ntsteps,t,'Parent', fig1); hold on;
pcolor(hfig1{t},X,Y,LStotal_temp); 
clim([0,1])
hold on;plot(hfig1{t},rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r')
hold on;plot(hfig1{t},rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
hold on;plot(hfig1{t},true_wpos(t,1),true_wpos(t,2),'r*','Markersize',10, 'LineWidth',3)
colorbar
xlabel(' x (m)'),ylabel('y (m)')
hold off

%//////////////////////////////////////////////////////////////////
% Apply Surface dilation filter (imdilate) to compensate for whale movement:
%dx=5; %Grid spacing SELCETED BASED ON WHAT??)
mean_horiz_swimspeed= 0.5; %mean horizontal swim speed (for sperm whale) [m/s]
dt90 = 9.73*abs(t-tstep90); %Elapsed time from when animals are at 90deg (0 tdoa) [s]
%(50m per time step equates to 9.73 s per time step if 
% the boat moves in each time step 50 m assuming it travels 10kts (5.144m/s))
if dt90==0
dt90=2;
end
mdwh=dt90*mean_horiz_swimspeed; %maximum distance whale could have travelled horizontally [m]

%gridpoints for filter in x (and also y since the same assumptions re swim speed and grid space)
xgridf=0:dx:(mdwh+dx);
filt_x_grid = [-fliplr(xgridf(2:end)),xgridf];
[Fx,Fy] = meshgrid(filt_x_grid,filt_x_grid);
De = Fx.^2/mdwh^2 + Fy.^2/mdwh^2; %normalize to max distance that the whale could have swam (feasible distance will be 1 or lower/ or is it 2 or lower??)
%F = zeros(size(De)); F((De <=1)) = 1;

LSdilate(:,:,t)=imdilate(LStotal_temp,(De<=1));

%//////////////////////////////////////////////////////////////////


end
%fontsize(fig1,14,'points')

% PLOT dilated surface (for some reason it messes up fig 1 if plotted
% within the loop)
fig2=figure; hfig2=cell(1,Ntsteps);
for t=1:Ntsteps
ip1=1;ip2=2;
rp=hyph_pos(:,:,t);
%PLOT
hfig2{t}=subplot(1,Ntsteps,t,'Parent', fig2); hold on;
pcolor(hfig2{t},X,Y,LSdilate(:,:,t)); hold on;
clim([0,1])
plot(hfig2{t},rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
plot(hfig2{t},rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
plot(hfig2{t},true_wpos(t,1),true_wpos(t,2),'r*','Markersize',10, 'LineWidth',3)
colorbar
xlabel(' x (m)'),ylabel('y (m)')
hold off
end


% //////////////////////// Normal LS without dilation//////////////////////
LStotal_temp=prod(LS,3);
LStotal=reshape(LStotal_temp,[Ngp_x,Ngp_x,Ngp_z]);
LStotal_hyperbolas_temp = sum(LS_Hyperbolas,3);
LStotal_hyperbolas = reshape(LStotal_hyperbolas_temp,[Ngp_x,Ngp_x,Ngp_z]);
fig3=figure;
pcolor(X,Y,LStotal),hold on
plot([hyph_pos(ip1,1,1),hyph_pos(ip1,1,end)],[hyph_pos(ip1,2,1),hyph_pos(ip1,2,end)],'r-', 'Linewidth', 3),hold on
for k=2:Ntsteps
plot([true_wpos(k-1,1),true_wpos(k,1)],[true_wpos(k-1,2),true_wpos(k,2)],'r*-','Markersize',10,'Linewidth', 1)
end
xlabel(' x (m)'),ylabel('y (m)')
clim([0,1])
colorbar
fontsize(fig3,14,'points')

fig4=figure;
imagesc(x,y,LStotal_hyperbolas); hold on
set(gca,'YDir', 'normal');
plot([hyph_pos(ip1,1,1),hyph_pos(ip1,1,end)],[hyph_pos(ip1,2,1),hyph_pos(ip1,2,end)],'r-', 'Linewidth', 3),hold on
for k=2:Ntsteps
plot([true_wpos(k-1,1),true_wpos(k,1)],[true_wpos(k-1,2),true_wpos(k,2)],'r*-','Markersize',10,'Linewidth', 1)
end
xlabel(' x (m)'),ylabel('y (m)')
clim([0,1])
colorbar
fontsize(fig4,14,'points')


% [LSval,indx]=max(LStotal, [], 'all');
% estimated_position= [X(indx),Y(indx),0]
% LSval
% true_end_position =true_wpos(end,:)

[peakdata] =findpeaks2D(X,Y,LStotal);
fprintf(['True whale end location is: [', repmat('%g, ', 1, numel(true_wpos(end,:))-1), '%g]\n'],true_wpos(end,:))
%take just one of the estimated peaks (the other is mirror image)
n=1;
estimated_position=[peakdata.peakX(n),peakdata.peakY(n),0]; %were at the moment considering 2D, so z coordinate = 0.
fprintf(['Estimated whale location is: [', repmat('%g, ', 1, numel(estimated_position)-1), '%g]\n'],estimated_position)
fprintf('Width of the peak in X dirextion %.2f \n',peakdata.peakXWidth(n))
fprintf('Width of the peak in Y dirextion %.2f \n',peakdata.peakYWidth(n))
fprintf('Ambiguity value for the estimated location is %.2f \n',peakdata.peakZ(n))
fprintf('The variance used for Gaussian is %.3f \n',sig)
fprintf('Number of sensors used is %.0f \n',Nsensors)
%//////////////////////////////////////////////////////////////////////////


%//////////////////////// LS with dilation//////////////////////
LStotaldilate=prod(LSdilate,3);
fig7=figure;
pcolor(X,Y,LStotaldilate),hold on
plot([hyph_pos(ip1,1,1),hyph_pos(ip1,1,end)],[hyph_pos(ip1,2,1),hyph_pos(ip1,2,end)],'r-', 'Linewidth', 3),hold on
for k=2:Ntsteps
plot([true_wpos(k-1,1),true_wpos(k,1)],[true_wpos(k-1,2),true_wpos(k,2)],'r*-','Markersize',10,'Linewidth', 1)
end
xlabel(' x (m)'),ylabel('y (m)')
clim([0,1])
colorbar
fontsize(fig7,14,'points')


[peakdata_dilate] =findpeaks2D(X,Y,LStotaldilate);
fprintf(['True whale end location is: [', repmat('%g, ', 1, numel(true_wpos(end,:))-1), '%g]\n'],true_wpos(end,:))
%take just one of the estimated peaks (the other is mirror image)
n=1;
estimated_position=[peakdata_dilate.peakX(n),peakdata_dilate.peakY(n),0]; %were at the moment considering 2D, so z coordinate = 0.
fprintf(['Estimated whale location is: [', repmat('%g, ', 1, numel(estimated_position)-1), '%g]\n'],estimated_position)
fprintf('Width of the peak in X dirextion %.2f \n',peakdata_dilate.peakXWidth(n))
fprintf('Width of the peak in Y dirextion %.2f \n',peakdata_dilate.peakYWidth(n))
fprintf('Ambiguity value for the estimated location is %.2f \n',peakdata_dilate.peakZ(n))
fprintf('The variance used for Gaussian is %.3f \n',sig)
fprintf('Number of sensors used is %.0f \n',Nsensors)

%//////////////////////////////////////////////////////////////////////////





%//////////////////////////////////////////////////////////////////////////
%% SIMULATION 4: Moving source, moving array with more sensors, no noise on measurements

clear,close all
%Determine sigma (std) for the Gaussian:
sig =0.003; %0.005 (if whale moves faster)

Ntsteps=4; %number of time steps

% Simulate Hydrophone positions - 3 sensors

hyph_pos(:,:,1)=[0,0,0;3,0,0;40,0,0]; %[x1,y1,z1; x2,y2,z2; x3,y3,z3];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[50,0,0;50,0,0;50,0,0];
end
Nsensors=size(hyph_pos,1);

% True whale position
true_wpos(1,:)= [100,50,0]; %[x,y,z];
for t=2:Ntsteps
true_wpos(t,:)= true_wpos(t-1,:) + [3,6,0]; %[x,y,z];  + [10,16,0]; %for faster moving whale (almost 4kts)
end

% True TDOAs
% compute tdoas between all sensors:
c=1500;
Ncombos=Nsensors*(Nsensors-1)- sum((Nsensors-1):-1:1); %Number of tdoa pairs
true_tdoas= cell(1,Ncombos);
count=1;
for ip1=1:Nsensors
    for ip2= ip1+1:Nsensors
        true_tdoa_temp=zeros(1,Ntsteps);
        for t=1:Ntsteps % at each time t
            rp=hyph_pos(:,:,t);

            dt1 = 1/c*sqrt((rp(ip1,1)-true_wpos(t,1))^2 + ...
                (rp(ip1,2)-true_wpos(t,2))^2 + ...
                (rp(ip1,3)-true_wpos(t,3))^2);
            dt2 = 1/c*sqrt((rp(ip2,1)-true_wpos(t,1))^2 + ...
                (rp(ip2,2)-true_wpos(t,2))^2 + ...
                (rp(ip2,3)-true_wpos(t,3))^2);
            true_tdoa_temp(t) = dt1-dt2;
        end
        true_tdoas{count}=true_tdoa_temp;
        count=count+1;
    end
end

%re-arrange so that the measured tdoas are per time step (matrix of (Ncombos,Ntsteps)):
tdoa_measured = vertcat(true_tdoas{:}); %Assume for the moment no noise in measurements

% Cell grids to evaluate:
% It could move with the hydrophone position or cover an entire area - what
% is better??
x=-200:5:200;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate

LStotal_tstep=zeros(1,N,Ntsteps);

for t=1:Ntsteps % for each time step compute LS (should be a product of LSs for all hydrophone pairs)

    rp=hyph_pos(:,:,t);
    count=1;

    LS_temp=zeros(Ncombos,N);
    figure, hold on
    for ip1=1:Nsensors
        for ip2= ip1+1:Nsensors
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

            LS_temp(count,:)= exp(-1/(2*sig^2).*(tdoa_model-tdoa_measured(count,t)).^2);
            

            LStotal_temp=reshape(LS_temp(count,:),[Ngp_x,Ngp_x,Ngp_z]);
            subplot(1,Nsensors,count)
            pcolor(X,Y,LStotal_temp); hold on
            clim([0,1])
            plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
            plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
            plot(true_wpos(t,1),true_wpos(t,2),'k*')
            colorbar
            xlabel(' x (m)'),ylabel('y (m)')
            title(['Ambiguity surface at time step ',num2str(t),...
                ', sensors ', num2str(ip1), ' and ', num2str(ip2)])

            count=count+1;

        end
    end

    LStotal_tstep(:,:,t) = prod(LS_temp);
    LStotal_tstep_temp=reshape(LStotal_tstep(:,:,t),[Ngp_x,Ngp_x,Ngp_z]);
    figure,
    pcolor(X,Y,LStotal_tstep_temp); hold on
    clim([0,1])
    for k=1:Nsensors
        plot(rp(k,1),rp(k,2),'r^','MarkerFaceColor','r'),hold on
    end
    plot(true_wpos(t,1),true_wpos(t,2),'k*')
    colorbar
    xlabel(' x (m)'),ylabel('y (m)')
    title(['Ambiguity surface at time step ',num2str(t)])

end

LStotal_temp=prod(LStotal_tstep,3);
LStotal=reshape(LStotal_temp,[Ngp_x,Ngp_x,Ngp_z]);
figure,
pcolor(X,Y,LStotal),hold on
hph=1;
plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-'),hold on
for k=2:Ntsteps
plot([true_wpos(k-1,1),true_wpos(k,1)],[true_wpos(k-1,2),true_wpos(k,2)],'k*-')
end
xlabel(' x (m)'),ylabel('y (m)'), title ('Total ambiguity surface')
clim([0,1])
colorbar

%this is same as below with findpeaks, but findpeaks reports the width of
%the peak as well
% [LSval,indx]=max(LStotal, [], 'all');
% estimated_position= [X(indx),Y(indx),0]
% LSval
% true_end_position =true_wpos(end,:)

findpeaks_2d_loop
fprintf(['True whale end location is: [', repmat('%g, ', 1, numel(true_wpos(end,:))-1), '%g]\n'],true_wpos(end,:))
%take just one of the estimated peaks (the other is mirror image)
n=1;
estimated_position=[peakdata.peakX(n),peakdata.peakY(n),0]; %were at the moment considering 2D, so z coordinate = 0.
fprintf(['Estimated whale location is: [', repmat('%g, ', 1, numel(estimated_position)-1), '%g]\n'],estimated_position)
fprintf('Width of the peak in X dirextion %.2f \n',peakdata.peakXWidth(n))
fprintf('Width of the peak in Y dirextion %.2f \n',peakdata.peakYWidth(n))
fprintf('Ambiguity value for the estimated location is %.2f \n',peakdata.peakZ(n))
fprintf('The variance used for Gaussian is %.3f \n',sig)
fprintf('Number of sensors used is %.0f \n',Nsensors)
% fprintf(' Estimated locations in X dirextion %.2f \n',peakdata.peakX)
% fprintf(' Estimated locations in Y dirextion %.2f \n',peakdata.peakY)



%//////////////////////////////////////////////////////////////////////////
%% SIMULATION 5: Moving array (2 sensors), two stationary sources, no noise

clear, close all
%Determine sigma (variance) for the Gaussian:
sig =0.003;

Ntsteps=4; %number of time steps

% Simulate Hydrophone positions
hyph_pos(:,:,1)=[0,0,0;40,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[50,0,0;50,0,0];
end

% True whales positions
true_wpos(1,:)= [100,50,0]; %[x,y,z];
true_wpos(2,:)= [10,100,0];
Nsources=size(true_wpos,1);

% True TDOAs
true_tdoa=zeros(Nsources,Ntsteps);
ip1=1;ip2=2;
c=1500;
for t=1:Ntsteps
    rp=hyph_pos(:,:,t);
    dt1 = 1/c*sqrt((rp(ip1,1)-true_wpos(:,1)).^2 + ...
        (rp(ip1,2)-true_wpos(:,2)).^2 + ...
        (rp(ip1,3)-true_wpos(:,3)).^2);
    dt2 = 1/c*sqrt((rp(ip2,1)-true_wpos(:,1)).^2 + ...
        (rp(ip2,2)-true_wpos(:,2)).^2 + ...
        (rp(ip2,3)-true_wpos(:,3)).^2);
    true_tdoa(:,t) = dt1-dt2;
end


tdoa_measured = true_tdoa; %Assume for the moment no noise in measurements
%measurements are a matrix Nsources x Nsteps

% Cell grids to evaluate:
% It could move with the hydrophone position or cover an entire area - what
% is better??
x=-200:5:200;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate



for t=1:Ntsteps % for each time step compute LS

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

LS(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-tdoa_measured(:,t)).^2);

figure,hold on
for n=1:Nsources
LStotal_temp=reshape(LS(n,:,t),[Ngp_x,Ngp_x,Ngp_z]);
subplot(1,Nsources,n)
pcolor(X,Y,LStotal_temp); hold on
clim([0,1])
plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
plot(true_wpos(n,1),true_wpos(n,2),'k*')
colorbar
xlabel(' x (m)'),ylabel('y (m)')
title(['LS for source ',num2str(n),' time step ', num2str(t)])
end

end



LStotal_temp=prod(LS,3);
figure, hold on
for n=1:Nsources
LStotal=reshape(LStotal_temp(n,:),[Ngp_x,Ngp_x,Ngp_z]);
subplot(1,Nsources,n)
pcolor(X,Y,LStotal),hold on
hph=1;
plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-'),hold on
plot(true_wpos(n,1),true_wpos(n,2),'r*')
xlabel(' x (m)'),ylabel('y (m)')
title ('Total ambiguity surface')
clim([0,1])
colorbar
end


%% REAL DATA example - Fragments of TDOA tracks
%Lasker AC 109 - processed with my TDOA tracking package (a chunk of
%encounter)
clear,close all
load('./LocalizationTest/LaskerAC109_Results/Lasker_AC109_Results.mat')

% Find all tracks that cross the beam - search within x degrees of the
% beam)
bear2tdoa = @(x) cosd(x).*(parameters.d/parameters.c);
maxtdoa=bear2tdoa(85); %or try 70 degrees to include more

indx=nan(size(Tracks,2),1);
for k=1:size(Tracks,2)
indx(k)=ismembertol(0,Tracks(k).tdoa,maxtdoa); %Checked and works good
end
indx=logical(indx);
Nsources=sum(indx);
Tracks_selected=Tracks(indx);

% Simulate array positions (assume the boat is traveling at 10kts=5.14 m/s)
mint=min(vertcat(Tracks_selected(:).time));
maxt=max(vertcat(Tracks_selected(:).time));
timevec=mint:parameters.dt:maxt;
Ntsteps=numel(timevec); %number of time steps

boat_move_m=5.14*parameters.dt;

% Simulate Hydrophone positions
hyph_pos(:,:,1)=[0,0,0;parameters.d,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[boat_move_m,0,0;boat_move_m,0,0];
end


% Re-arrange tracks in a matrix (Nsources x Ntsteps)
tdoa_measured= nan(Nsources,Ntsteps);
for k=1:Nsources
    time_indx=ismember(timevec,Tracks_selected(k).time);
    tdoa_measured(k,time_indx)= Tracks_selected(k).tdoa;
end


% Compute modeled TDOAs
% Create a grid
x=0:10:3500;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate
% Specify sound speed
c=1500;
% Specify variance for the Gaussian (should reflect all errors):
% sig=0.0024; % From Yvonnes paper
sig=0.008;


LS= nan(Nsources,N,Ntsteps);

for t=1:Ntsteps % for each time step compute LS

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

LS(:,:,t)= exp(-1/(2*sig^2).*(tdoa_model-(-1.*tdoa_measured(:,t))).^2);

% figure,hold on
% for n=1:Nsources
% LStotal_temp=reshape(LS(n,:,t),[Ngp_x,Ngp_x,Ngp_z]);
% subplot(1,Nsources,n)
% pcolor(X,Y,LStotal_temp); hold on
% clim([0,1])
% plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
% plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
% plot(true_wpos(n,1),true_wpos(n,2),'k*')
% colorbar
% xlabel(' x (m)'),ylabel('y (m)')
% title(['LS for source ',num2str(n),' time step ', num2str(t)])
% end

end


LStotal_temp=prod(LS,3,'omitnan');
for n=1:Nsources
    figure, hold on
    LStotal=reshape(LStotal_temp(n,:),[Ngp_x,Ngp_x,Ngp_z]);
    s=pcolor(X,Y,LStotal);
    s.EdgeColor='none';
    hold on
    hph=1;
    plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-'),hold on
    xlabel(' x (m)'),ylabel('y (m)')
    title ('Total ambiguity surface for source ', num2str(n))
    clim([0,1])
    colorbar
end

%% From Track 1 I get a small AS value - probably animals move too much in
%that time- select narrower time when they are around the beam (and use
%bigger variance than what is used in Yvonnes paper)
clear,close all
addpath('/Users/pinagruden/Dropbox/Pina/Git_projects/TDOA_tracking_master'); %for plotting

load('./LocalizationTest/LaskerAC109_Results/Lasker_AC109_Results.mat')
load('./LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_clicks_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL')
Rxy_envelope_both{1}=Rxy_envelope_ALL.*scalar_clicks;
load('./LocalizationTest/LaskerAC109_Raw_CrossCorrelogram/Lasker_AC109_whistles_rawCrossCorrelogram_ALL.mat',...
    'Rxy_envelope_ALL','lags','t_serialdate')
Rxy_envelope_both{2}=Rxy_envelope_ALL.*scalar_whistles;


% ////// Find all tracks that cross the beam /////////////

%1. Search for all tracks within x degrees of the beam
bear2tdoa = @(x) cosd(x).*(parameters.d/parameters.c);
maxtdoa=bear2tdoa(60);

indx=nan(size(Tracks,2),1);
for k=1:size(Tracks,2)
indx(k)=ismembertol(0,Tracks(k).tdoa,maxtdoa); %Checked and works good
end
indx=logical(indx);
Nsources=sum(indx);
Tracks_selected=Tracks(indx);

%2. Identify min, max times for these tracks & create a time vector
mint=min(vertcat(Tracks_selected(:).time));
maxt=max(vertcat(Tracks_selected(:).time));
timevec=mint:parameters.dt:maxt;
Ntsteps=numel(timevec); %number of time steps


%3. PLOT All tracks and selected:
parameters.saveplotsofresults=0;
measure.T=0;
plot_results(t_serialdate,lags,Rxy_envelope_both,measure,Tracks, parameters)
figure(2), hold on
for k=1:Nsources
text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1),num2str(k),'Color','k')
end



% ////////// Simulate Boat movement and Hydrophone positions ///////////
%2. Simulate Hydrophone positions

boat_move_m=5.14*parameters.dt; %(assume the boat is traveling at 10kts=5.14 m/s)

hyph_pos(:,:,1)=[0,0,0;parameters.d,0,0]; %[x1,y1,z1; x2,y2,z2];
for t=2:Ntsteps
hyph_pos(:,:,t)=hyph_pos(:,:,t-1)+[boat_move_m,0,0;boat_move_m,0,0];
end


%/////////// Re-arrange tracks in a matrix (Nsources x Ntsteps) ///////////
tdoa_measured= nan(Nsources,Ntsteps);
for k=1:Nsources
    time_indx=ismember(timevec,Tracks_selected(k).time);
    tdoa_measured(k,time_indx)= Tracks_selected(k).tdoa;
end


%////// Compute modeled TDOAs//////
%1. Create a grid
x=0:10:4000;
Ngp_x=length(x);
y=x;
z=0;
Ngp_z=length(z);
[X,Y,Z] = meshgrid(x,y,z);
wpos=[X(:),Y(:),Z(:)]; %gives grid of N x [x,y,z] coordinates
N= size(wpos,1); %number of grid points to evaluate

%2. Specify sound speed
c=1500;

%3. Specify variance for the Gaussian (should reflect all errors):
% sig=0.0024; % From Yvonnes paper
sig=0.008;


%% ///////////////  Compute LS for a selected source ////////////////////
% See the plotted Figure 2 and Select track/s you want to look at:

% % Select one track:
% SelectedTrack =13;
% ind=find(abs(Tracks_selected(SelectedTrack).tdoa)<maxtdoa);
% selected_indx=ismember(timevec,Tracks_selected(SelectedTrack).time(ind));
% tdoa_measured_select=tdoa_measured(SelectedTrack,:);

%Or combine tracks:
SelectedTracks =[13,16,20,23];
selected_indx=zeros(length(SelectedTracks),length(timevec));
for n=1:length(SelectedTracks)
ind=find(abs(Tracks_selected(SelectedTracks(n)).tdoa)<maxtdoa);
selected_indx(n,:)=ismember(timevec,Tracks_selected(SelectedTracks(n)).time(ind));
end
selected_indx=sum(selected_indx,1);
tdoa_measured_select=sum(tdoa_measured(SelectedTracks,:),1,'omitnan');
tdoa_measured_select(all(isnan(tdoa_measured(SelectedTracks,:)),1)) = NaN;



% Pre-allocate
LS_select= nan(1,N,Ntsteps);

count=1; plotf=0;%[1,50,100,150,200,250];
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
            LStotal_temp=reshape(LS_select(:,:,t),[Ngp_x,Ngp_x,Ngp_z]);
            s=pcolor(X,Y,LStotal_temp);
            s.EdgeColor='none';
            hold on
            clim([0,1])
            plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
            plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
            colorbar
            xlabel(' x (m)'),ylabel('y (m)')
            title(['LS for source 1, time step ', num2str(t)])
        end

        count=count+1;
    end
end

LStotal_temp=prod(LS_select,3,'omitnan');
figure, hold on
LStotal=reshape(LStotal_temp(1,:),[Ngp_x,Ngp_x,Ngp_z]);
s=pcolor(X,Y,LStotal);
s.EdgeColor='none';
hold on
hph=1;
plot([hyph_pos(hph,1,1),hyph_pos(hph,1,end)],[hyph_pos(hph,2,1),hyph_pos(hph,2,end)],'r-'),hold on
xlabel(' x (m)'),ylabel('y (m)')
title (['Total ambiguity surface for source ', num2str(SelectedTracks)])
clim([0,1])
colorbar
xlim([0,4000])
set(gca,'FontSize',16)


%% Connect fragments


figure(2), hold on
for k=1:Nsources
text(Tracks_selected(k).time_local(1),Tracks_selected(k).tdoa(1),num2str(k),'Color','k')
end

coeffs=zeros(Nsources,2);
for k=1:Nsources
coeffs(k,:)=polyfit(Tracks_selected(k).time,Tracks_selected(k).tdoa,1);
end

figure,
plot(coeffs(:,1),coeffs(:,2),'k*'), xlabel('slope'),ylabel('intercept')
for k=1:Nsources
text(coeffs(k,1),coeffs(k,2),num2str(k),'Color','r')
end

%%

figure(2), hold on
for k=1:M
text(Tracks(k).time_local(1),Tracks(k).tdoa(1),num2str(k),'Color','b')
end

M=size(Tracks,2);
coeffs=zeros(M,2);
for k=1:M
coeffs(k,:)=polyfit(Tracks(k).time,Tracks(k).tdoa,1);
end

figure,
plot(coeffs(:,1),coeffs(:,2),'k*'), xlabel('slope'),ylabel('intercept')
for k=1:M
text(coeffs(k,1),coeffs(k,2),num2str(k),'Color','r')
end