% Localize with Ambiguity surfaces


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Make sure you had first run A0_SetUp.m

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



%% /////////////////////1) Select TDOA track///////////////////
% Select which TDOA track or which fragments you want to compute
% localization for

[tdoa_measured_select, selected_indx] = select_tracks(Tracks_selected,timevec,AS_params.tdoa_cutoff);

%% /////////////////////Compute Ambiguity Surfaces///////////////////
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