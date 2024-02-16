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


Ntsteps = size(tdoa_measured_select,2); %number of time steps

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


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2.b) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
[AStotal_rough,ASdilatetotal_rough] = obtain_ambiguitysurface(hyph_pos,selected_indx, ...
    AS_params,BA_params,gridparams,wpos2D,tdoa_measured_select,sig,tstep90);


%~~~~~~~2.d) Determine the peak (Loc) in AS to get narrower x/yrange~~~~~~~
%Compute it for both normal AS and dilated and take whathever is smaller&bigger

%----------------------------------------------------------------
%----------------------------------------------------------------
% Select tallest peak as localization (since it's biambigous)
% ----------------Non-dilated surfaces-------------------
% Marginalize along x-axis
x_marg= sum(AStotal_rough,1);
[height_x]=findpeaks(x_marg,gridparams.X(1,:),'SortStr', 'descend', 'NPeaks', 1); 
% Take 50% from the peak to be my bounds
ind_xmarg=find(x_marg>height_x*0.5);

% If peak is extremely narrow then sometimes not both sides from the peak
% are included- ensure there are at least 3 values in there with the peak
% value in the middle:
if numel(ind_xmarg)<3
    ind_xpeak=find(x_marg==height_x);
    ind_xmarg = [ind_xpeak-1,ind_xpeak,ind_xpeak+1];
end

% Marginalize along y-axis
y_marg=sum(AStotal_rough,2);
[height_y]=findpeaks(y_marg,gridparams.Y(:,1),'SortStr', 'descend', 'NPeaks', 1);

% If there are insufficient bearing changess to have
% a peak in the surface- it will not have a peak so just take the max:
if isempty(height_y)
    height_y=max(y_marg);
    % Take 80% from the peak to be my bounds
    ind_ymarg=find(y_marg>height_y*0.2);
else % there is a defined peak in y marginalized surface
    % Take 50% from the peak to be my bounds
    ind_ymarg=find(y_marg>height_y*0.5);
end

% If peak is extremely narrow then sometimes not both sides from the peak
% are included- ensure there are at least 3 values in there with the peak
% value in the middle (although very unlikely in the y-direction)
if numel(ind_ymarg)<3
    ind_ypeak=find(y_marg==height_y);
    ind_ymarg = [ind_ypeak-1,ind_ypeak,ind_ypeak+1];
end

xrange_new_temp1=round([gridparams.x(ind_xmarg(1))-gridparams.dx,gridparams.x(ind_xmarg(end))+gridparams.dx]);
yrange_new_temp1=round([gridparams.y(ind_ymarg(1))-gridparams.dy,gridparams.y(ind_ymarg(end))+gridparams.dy]);

%------------------ Dilated surfaces----------------------
% Marginalize along x-axis
x_marg_dilate= sum(ASdilatetotal_rough,1);
[height_x]=findpeaks(x_marg_dilate,gridparams.X(1,:),'SortStr', 'descend', 'NPeaks', 1);
% Take 50% from the peak to be my bounds
ind_xmarg=find(x_marg_dilate>height_x*0.5);

% If peak is extremely narrow then sometimes not both sides from the peak
% are included- ensure there are at least 3 values in there with the peak
% value in the middle:
if numel(ind_xmarg)<3
    ind_xpeak=find(x_marg_dilate==height_x);
    ind_xmarg = [ind_xpeak-1,ind_xpeak,ind_xpeak+1];
end

% Marginalize along y-axis
y_marg_dilate=sum(ASdilatetotal_rough,2);
[height_y]=findpeaks(y_marg_dilate,gridparams.Y(:,1), 'SortStr', 'descend', 'NPeaks', 1);

% If there are insufficient bearing changess to have
% a peak in the surface- it will not have a peak so just take the max:
if isempty(height_y)
    height_y=max(y_marg_dilate);
    % Take 80% from the peak to be my bounds
    ind_ymarg=find(y_marg_dilate>height_y*0.2);
else % there is a defined peak in y marginalized surface
    % Take 50% from the peak to be my bounds
    ind_ymarg=find(y_marg_dilate>height_y*0.5);
end

% If peak is extremely narrow then sometimes not both sides from the peak
% are included- ensure there are at least 3 values in there with the peak
% value in the middle (although very unlikely in the y-direction)
if numel(ind_ymarg)<3
    ind_ypeak=find(y_marg_dilate==height_y);
    ind_ymarg = [ind_ypeak-1,ind_ypeak,ind_ypeak+1];
end

xrange_new_temp2=round([gridparams.x(ind_xmarg(1))-gridparams.dx,gridparams.x(ind_xmarg(end))+gridparams.dx]);
yrange_new_temp2=round([gridparams.y(ind_ymarg(1))-gridparams.dy,gridparams.y(ind_ymarg(end))+gridparams.dy]);

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3.b) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
[AStotal,ASdilatetotal,AStotal_hyperbolas] = obtain_ambiguitysurface(hyph_pos,selected_indx, ...
    AS_params,BA_params,gridparams,wpos2D,tdoa_measured_select,sig,tstep90);

%add total grid- wpos2D to the gridparams:
gridparams.wpos = wpos2D;


end