function [AStotal,ASdilatetotal,AStotal_hyperbolas,gridparams] = computeAS(tdoa_measured_select,selected_indx,hyph_pos,AS_params,BA_params)
%computeAS.m computes an ambiguity surface (AS) for a given source based on 
%the measured and modeled time difference of arrivals (TDOAs)
%
% INPUTS:
% - tdoa_measured_select : 1 x T vector containing tdoas of the selected source, 
%                 where T is number of time steps
% - selected_indx : 1 x T vector indicating in which time step target
%                   exists (1) and in which it does not (0), where T is 
%                   number of time steps
% - hyph_pos : a NxMxT array, where T is number of time steps, N=number of
%              sensors and M is number of coordinates (e.g. x,y or x,y,z)
% - AS_params : a structure containing parameters for ambiguity surface
%               computation
% - BA_params : a structure containing parameters for boat and array
%
% OUTPUTS:
% - AStotal : final ambiguity surface, an A x B matrix where A is number of
%              y coordinates, and B is number of x coordinates 
% - ASdilatetotal: final ambiguity surface with dilation, same dimensions
%                  as AStotal
% - AStotal_hyperbolas: final surface with intersecting hyperbolas, same
%                       dimensions as AStotal
%
%
% Pina Gruden, Dec 2022, UH Manoa


% wpos=AS_params.wpos;
% Ngp_x=AS_params.Ngp_x;
% Ngp_y=AS_params.Ngp_y;
% dx=AS_params.dx;
% dy=AS_params.dy;
% sig=AS_params.sig;
sig_hyperbolas=AS_params.sig_hyperbolas;
c= AS_params.c;
mean_horiz_swimspeed = AS_params.mean_horiz_swimspeed;
timestep = BA_params.timestep;

Ntsteps = size(tdoa_measured_select,2); %number of time steps
% N= size(wpos,1); %number of grid points to evaluate

%~~~~~~~ 1) Get time step when animal is closest to 90deg (tdoa=0s)~~~~~~~~
% This is for the purpose of Surface Dilation
[~,ind]=min(abs(tdoa_measured_select));
tsteps=1:1:Ntsteps;
tstep90=tsteps(ind); % time step whenanimal is closest to the beam 

%////////////////// 2) FIRST COMPUTE AS with ROUGH GRID /////////////////

%~~~~~~~~~~~~~~~~~~~~~~~~~ 2.a) Create ROUGH GRID ~~~~~~~~~~~~~~~~~~~~~~~~~
dx=100; %in meters- make it very coarse
dy=100;
sig=AS_params.sig*10; % needs to be bigger since we're using coarser resolution
xrange=AS_params.xrange;
yrange=AS_params.yrange;
[gridparams,wpos2D] = make_2Dgrids(xrange,yrange,dx,dy);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2.b) Pre-allocate ~~~~~~~~~~~~~~~~~~~~~~~~~~
AS_select_rough= nan(gridparams.Ngp_y,gridparams.Ngp_x,Ntsteps);%swap x and y inpuput arguments since 
% reshape() will be used to fill in the values
ASdilate_rough = AS_select_rough;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2.c) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
%count=1; plotf=0;
for t=1:Ntsteps % for each time step compute LS
    if selected_indx(t)
        rp=hyph_pos(:,:,t);
        ip1=1;ip2=2;
        % for each hypothetical grid position wpos(wpi) this calculates
        % time between position and hydrophone position - dt1 (between
        % hyph 1 and position wpos(wpi)), dt2 (between hyph2 and position wpos(wpi)), then
        % tdoa_model is tdoa between the hyph1 and hyph2 if the source
        % is at wpos(wpi).
        dt1= 1/c.*sqrt((rp(ip1,1)-wpos2D(:,1)).^2 +(rp(ip1,2)-wpos2D(:,2)).^2);
        dt2= 1/c.*sqrt((rp(ip2,1)-wpos2D(:,1)).^2 +(rp(ip2,2)-wpos2D(:,2)).^2);
        tdoa_model=dt1-dt2;
        tdoa_model=reshape(tdoa_model,[gridparams.Ngp_y,gridparams.Ngp_x]);

        tdoa_diff=(tdoa_model-tdoa_measured_select(:,t)).^2;
        AS_select_rough(:,:,t)= exp(-1/(2*sig^2).*tdoa_diff);
     
        %///////////////////DILATE SURFACES////////////////////////////////
        % Apply Surface dilation filter (imdilate) to compensate for whale movement:
        AStotal_temp=AS_select_rough(:,:,t);
        dt90 = timestep*abs(t-tstep90); %Elapsed time from when animals are at 90deg (0 tdoa) [s]
        if dt90==0

            ASdilate_rough(:,:,t)=AStotal_temp; %do not dilate if animal is at the beam (the tdoas should be correct)

        else
            mdwh=dt90*mean_horiz_swimspeed; %maximum distance whale could have travelled horizontally [m]

            %gridpoints for filter in x (and also y since the same assumptions re swim speed and grid space)
            xgridf=0:dx:(mdwh+dx);
            filt_x_grid = [-fliplr(xgridf(2:end)),xgridf];
            ygridf=0:dy:(mdwh+dy);
            filt_y_grid = [-fliplr(ygridf(2:end)),ygridf];
            [Fx,Fy] = meshgrid(filt_x_grid,filt_y_grid);
            De = sqrt(Fx.^2+Fy.^2)/mdwh; %normalize to max distance that the whale could have swam (feasible distance will be 1 or lower)

            ASdilate_rough(:,:,t)=imdilate(AStotal_temp,(De<=1));
        end
        %//////////////////////////////////////////////////////////////////

    end
end
clear AStotal_temp

% Get total AS without dilation:
AStotal_rough=prod(AS_select_rough,3,'omitnan');

% Get total AS with dilation:
ASdilatetotal_rough=prod(ASdilate_rough,3,'omitnan');

%~~~~~~~2.d) Determine the peak (Loc) in AS to get narrower x/yrange~~~~~~~
%Compute it for both normal AS and dilated and take whathever is smaller&bigger

%----------------------------------------------------------------
%----------------------------------------------------------------
% Select tallest peak as localization (since it's biambigous)
% ----------------Non-dilated surfaces-------------------
% Marginalize along x-axis
x_marg= sum(AStotal_rough,1);
[height_x]=findpeaks(x_marg,gridparams.X(1,:),'SortStr', 'descend', 'NPeaks', 1); 

% Marginalize along y-axis
y_marg=sum(AStotal_rough,2);
[height_y]=findpeaks(y_marg,gridparams.Y(:,1),'SortStr', 'descend', 'NPeaks', 1);

% If there are insufficient bearing changess to have
% a peak in the surface- it will not have a peak so just take the max:
if isempty(height_y) 
    height_y=max(y_marg);
end

% Take 50% from the peak to be my bounds
ind_xmarg=find(x_marg>height_x*0.5);
ind_ymarg=find(y_marg>height_y*0.5);

xrange_new_temp1=round([gridparams.x(ind_xmarg(1)),gridparams.x(ind_xmarg(end))]);
yrange_new_temp1=round([gridparams.y(ind_ymarg(1)),gridparams.y(ind_ymarg(end))]);

%------------------ Dilated surfaces----------------------
% Marginalize along x-axis
x_marg_dilate= sum(ASdilatetotal_rough,1);
[height_x]=findpeaks(x_marg_dilate,gridparams.X(1,:),'SortStr', 'descend', 'NPeaks', 1);

% Marginalize along y-axis
y_marg_dilate=sum(ASdilatetotal_rough,2);
[height_y]=findpeaks(y_marg_dilate,gridparams.Y(:,1), 'SortStr', 'descend', 'NPeaks', 1);

% If there are insufficient bearing changess to have
% a peak in the surface- it will not have a peak so just take the max:
if isempty(height_y) 
    height_y=max(y_marg);
end


% Take 50% from the peak to be my bounds
ind_xmarg=find(x_marg_dilate>height_x*0.5);
ind_ymarg=find(y_marg_dilate>height_y*0.5);

xrange_new_temp2=round([gridparams.x(ind_xmarg(1)),gridparams.x(ind_xmarg(end))]);
yrange_new_temp2=round([gridparams.y(ind_ymarg(1)),gridparams.y(ind_ymarg(end))]);

%-------------- Get the new range --------------------
xrange_new=[min(xrange_new_temp1(1),xrange_new_temp2(1)),max(xrange_new_temp1(2),xrange_new_temp2(2))];
yrange_new=[min(yrange_new_temp1(1),yrange_new_temp2(1)),max(yrange_new_temp1(2),yrange_new_temp2(2))];

%----------------------------------------------------------------
%----------------------------------------------------------------


%//////////////// 3) COMPUTE AS with FINER (desired) GRID /////////////////

%~~~~~~~~~~~~~~~~~~~~~~~~ 3.a) Create FINE GRID ~~~~~~~~~~~~~~~~~~~~~~~~~~~
dx=AS_params.dx; %now take the user specified resolution
dy=AS_params.dy;
sig=AS_params.sig; % user specified
xrange=xrange_new;
yrange=yrange_new;
[gridparams,wpos2D] = make_2Dgrids(xrange,yrange,dx,dy);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3.b) Pre-allocate ~~~~~~~~~~~~~~~~~~~~~~~~~~
AS_select= nan(gridparams.Ngp_y,gridparams.Ngp_x,Ntsteps);%swap x and y inpuput arguments since 
% reshape() will be used to fill in the values
AS_Hyperbolas = AS_select;
ASdilate = AS_select;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2.c) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
for t=1:Ntsteps % for each time step compute LS
    if selected_indx(t)
        rp=hyph_pos(:,:,t);
        ip1=1;ip2=2;
        % for each hypothetical grid position wpos(wpi) this calculates
        % time between position and hydrophone position - dt1 (between
        % hyph 1 and position wpos(wpi)), dt2 (between hyph2 and position wpos(wpi)), then
        % tdoa_model is tdoa between the hyph1 and hyph2 if the source
        % is at wpos(wpi).
        dt1= 1/c.*sqrt((rp(ip1,1)-wpos2D(:,1)).^2 +(rp(ip1,2)-wpos2D(:,2)).^2);
        dt2= 1/c.*sqrt((rp(ip2,1)-wpos2D(:,1)).^2 +(rp(ip2,2)-wpos2D(:,2)).^2);
        tdoa_model=dt1-dt2;
        tdoa_model=reshape(tdoa_model,[gridparams.Ngp_y,gridparams.Ngp_x]);

        tdoa_diff=(tdoa_model-tdoa_measured_select(:,t)).^2;
        AS_select(:,:,t)= exp(-1/(2*sig^2).*tdoa_diff);
        AS_Hyperbolas(:,:,t)= exp(-1/(2*sig_hyperbolas^2).*tdoa_diff);

        %///////////////////DILATE SURFACES////////////////////////////////
        % Apply Surface dilation filter (imdilate) to compensate for whale movement:
        AStotal_temp=AS_select(:,:,t);
        dt90 = timestep*abs(t-tstep90); %Elapsed time from when animals are at 90deg (0 tdoa) [s]
        if dt90==0

            ASdilate(:,:,t)=AStotal_temp; %do not dilate if animal is at the beam (the tdoas should be correct)

        else
            mdwh=dt90*mean_horiz_swimspeed; %maximum distance whale could have travelled horizontally [m]

            %gridpoints for filter in x (and also y since the same assumptions re swim speed and grid space)
            xgridf=0:dx:(mdwh+dx);
            filt_x_grid = [-fliplr(xgridf(2:end)),xgridf];
            ygridf=0:dy:(mdwh+dy);
            filt_y_grid = [-fliplr(ygridf(2:end)),ygridf];
            [Fx,Fy] = meshgrid(filt_x_grid,filt_y_grid);
            De = sqrt(Fx.^2+Fy.^2)/mdwh; %normalize to max distance that the whale could have swam (feasible distance will be 1 or lower)

            ASdilate(:,:,t)=imdilate(AStotal_temp,(De<=1));
        end
        %//////////////////////////////////////////////////////////////////

    end
end
clear AStotal_temp

% Get total AS without dilation:
AStotal=prod(AS_select,3,'omitnan');

% Get total AS with dilation:
ASdilatetotal=prod(ASdilate,3,'omitnan');

% Get intersecting hyperbolas:
AStotal_hyperbolas = sum(AS_Hyperbolas,3,'omitnan');

%add total grid- wpos2D to the gridparams:
gridparams.wpos = wpos2D;


end