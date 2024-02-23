function [AStotal,ASdilatetotal,AStotal_hyperbolas,gridparams,flag_locnotpossible] = computeAS(tdoa_measured_select,selected_indx,hyph_pos,AS_params,BA_params)
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
sig=AS_params.sig*3; % needs to be bigger since we're using coarser resolution
xrange=AS_params.xrange;
yrange=AS_params.yrange;
[gridparams,wpos2D] = make_2Dgrids(xrange,yrange,dx,dy);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2.b) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
[AStotal_rough,ASdilatetotal_rough] = obtain_ambiguitysurface(hyph_pos,selected_indx, ...
    AS_params,BA_params,gridparams,wpos2D,tdoa_measured_select,sig,tstep90);


%~~~~~~~2.c) Determine the peak (Loc) in AS to get narrower x/yrange~~~~~~~
%Compute it for both normal AS and dilated and take whathever is smaller&bigger

%----------------------------------------------------------------
% ----------------Non-dilated surfaces-------------------
[xrange_new_temp1,yrange_new_temp1,flag_locnotpossible1] = get_xyrange(AStotal_rough,gridparams);

%------------------ Dilated surfaces----------------------
[xrange_new_temp2,yrange_new_temp2,flag_locnotpossible2] = get_xyrange(ASdilatetotal_rough,gridparams);

if flag_locnotpossible1 || flag_locnotpossible2 %There is not enough change in bearing to get localization
    
    flag_locnotpossible = true;

    AStotal=[];
    ASdilatetotal=[];
    AStotal_hyperbolas=[];
    gridparams=[];

else %There is enough change in bearing to get localization

    flag_locnotpossible = false;

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

end