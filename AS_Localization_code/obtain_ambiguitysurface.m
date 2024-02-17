function [AStotal,ASdilatetotal,AStotal_hyperbolas] = obtain_ambiguitysurface(hyph_pos,selected_indx, ...
    AS_params,BA_params,gridparams,wpos2D,tdoa_measured_select,sig,tstep90)
% function obtain_ambiguitysurface.m computes an ambiguity surface for a
% given source for a specified grid and surface paramters.
%
% INPUTS:
% - hyph_pos : a NxMxT array, where T is number of time steps, N=number of
%              sensors and M is number of coordinates (e.g. x,y or x,y,z)
% - selected_indx : 1 x T vector indicating in which time step target
%                   exists (1) and in which it does not (0), where T is 
%                   number of time steps
% - AS_params : a structure containing parameters for ambiguity surface
%               computation
% - BA_params : a structure containing parameters for boat and array
% - gridparams : a structure contining grid parameters, such as x&y-ranges,
%               resolution, number of grid points
% - wpos2D : a N x 2 matrix, where N is number of all grid points and the
%           columns correspond to x and y coordinates respectively
% - tdoa_measured_select : 1 x T vector containing tdoas of the selected source, 
%                 where T is number of time steps
% - sig : sigma (standard deviation) for the Gaussian (used to create
%         ambiguity surface)
% - tstep90 : time step when animal is closest to 90deg
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
% Pina Gruden, Feb 2024, UH Manoa

c= AS_params.c;
mean_horiz_swimspeed = AS_params.mean_horiz_swimspeed;
sig_hyperbolas=AS_params.sig_hyperbolas;
timestep = BA_params.timestep;
Ntsteps = size(tdoa_measured_select,2); %number of time steps
dx=gridparams.dx;
dy=gridparams.dy;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Pre-allocate ~~~~~~~~~~~~~~~~~~~~~~~~~~
AS_select= nan(gridparams.Ngp_y,gridparams.Ngp_x,Ntsteps);%swap x and y inpuput arguments since 
% reshape() will be used to fill in the values
ASdilate = AS_select;
AS_Hyperbolas = AS_select;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
 
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

% Get total AS without dilation:
AStotal=prod(AS_select,3,'omitnan');

% Get total AS with dilation:
ASdilatetotal=prod(ASdilate,3,'omitnan');

% Get intersecting hyperbolas:
AStotal_hyperbolas = sum(AS_Hyperbolas,3,'omitnan');

end