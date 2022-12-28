function [AStotal,ASdilatetotal,AStotal_hyperbolas] = computeAS(tdoa_measured_select,selected_indx,hyph_pos,AS_params,BA_params)
%computeAS.m computes an ambiguity surface (AS) for a given source based on 
%the measured and modeled time difference of arrivals (TDOAs)

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


% OUTPUTS:
% - AStotal : final ambiguity surface, an A x B matrix where A is number of
%              y coordinates, and B is number of x coordinates 
% - ASdilatetotal: final ambiguity surface with dilation, same dimensions
%                  as AStotal
% - AStotal_hyperbolas: final surface with intersecting hyperbolas, same
%                       dimensions as AStotal

% Pina Gruden, Dec 2022, UH Manoa


wpos=AS_params.wpos;
Ngp_x=AS_params.Ngp_x;
Ngp_y=AS_params.Ngp_y;
dx=AS_params.dx;
dy=AS_params.dy;
sig=AS_params.sig;
sig_hyperbolas=AS_params.sig_hyperbolas;
c= AS_params.c;
mean_horiz_swimspeed = AS_params.mean_horiz_swimspeed;
timestep = BA_params.timestep;

Ntsteps = size(tdoa_measured_select,2); %number of time steps
N= size(wpos,1); %number of grid points to evaluate

%~~~~~~~ 1) Get time step when animal is closest to 90deg (tdoa=0s)~~~~~~~~
% This is for the purpose of Surface Dilation
[~,ind]=min(abs(tdoa_measured_select));
tsteps=1:1:Ntsteps;
tstep90=tsteps(ind); % time step whenanimal is closest to the beam 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2) Pre-allocate ~~~~~~~~~~~~~~~~~~~~~~~~~~
AS_select= nan(1,N,Ntsteps);
AS_Hyperbolas = AS_select;
ASdilate = nan(Ngp_y,Ngp_x,Ntsteps); %swap x and y inpuput arguments since 
% reshape() will be used to fill in the values


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3) Compute AS ~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        dt1= 1/c.*sqrt((rp(ip1,1)-wpos(:,1)).^2 +(rp(ip1,2)-wpos(:,2)).^2);
        dt2= 1/c.*sqrt((rp(ip2,1)-wpos(:,1)).^2 +(rp(ip2,2)-wpos(:,2)).^2);
        tdoa_model=dt1-dt2;

        tdoa_diff=(tdoa_model-tdoa_measured_select(:,t)).^2;
        AS_select(:,:,t)= exp(-1/(2*sig^2).*tdoa_diff);
        AS_Hyperbolas(:,:,t)= exp(-1/(2*sig_hyperbolas^2).*tdoa_diff);

        %///////////////////DILATE SURFACES////////////////////////////////
        % Apply Surface dilation filter (imdilate) to compensate for whale movement:
        AStotal_temp=reshape(AS_select(:,:,t),[Ngp_y,Ngp_x]);
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

%         %Plot surfaces:
%         if any(count==plotf)
%             figure,hold on
%             LStotal_temp=reshape(LS_select(:,:,t),[Ngp_y,Ngp_x]);
%             s=pcolor(X,Y,LStotal_temp);
%             s.EdgeColor='none';
%             clim([0,1])
%             axis equal
%             plot(rp(ip1,1),rp(ip1,2),'r^','MarkerFaceColor','r'),hold on
%             plot(rp(ip2,1),rp(ip2,2),'r^','MarkerFaceColor','r')
%             colorbar
%             xlabel(' x (m)'),ylabel('y (m)')
%             title(['LS for selected source, time step ', num2str(t)])
%         end
% 
%         count=count+1;

    end
end
clear AStotal_temp

% Get total AS without dilation:
AStotal_temp=prod(AS_select,3,'omitnan');
AStotal=reshape(AStotal_temp,[Ngp_y,Ngp_x]);

% Get total AS with dilation:
ASdilatetotal=prod(ASdilate,3,'omitnan');

% Get intersecting hyperbolas:
AStotal_hyperbolas_temp = sum(AS_Hyperbolas,3,'omitnan');
AStotal_hyperbolas = reshape(AStotal_hyperbolas_temp,[Ngp_y,Ngp_x]);

end