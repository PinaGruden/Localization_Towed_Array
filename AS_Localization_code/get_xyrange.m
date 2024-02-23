function [xrange_new,yrange_new,flag_locnotpossible] = get_xyrange(AStotal,gridparams)
% INPUTS:
% - AStotal : ambiguity surface, an A x B matrix where A is number of
%              y coordinates, and B is number of x coordinates  
% - gridparams : a structure contining grid parameters, such as x&y-ranges,
%               resolution, number of grid points
%
% OUTPUTS:
% -xrange_new - 1 x 2 vector of range [min, max] for x-coordinate
% -yrange_new - 1 x 2 vector of range [min, max] for y-coordinate
% -flag_locnotpossible - logical (true or false) if localization is
%       possible or not
%
% Pina Gruden, Feb 20224, UH Manoa

 
% Marginalize along x-axis
% x_marg= sum(AStotal,1);
x_marg= max(AStotal,[],1);
[height_x]=findpeaks(x_marg,gridparams.X(1,:),'SortStr', 'descend', 'NPeaks', 1); 

% Marginalize along y-axis
% y_marg=sum(AStotal,2);
y_marg=max(AStotal,[],2); %note taking max here since marginalization 
% (summing) somethimes does not produce a peak
[height_y]=findpeaks(y_marg,gridparams.Y(:,1),'SortStr', 'descend', 'NPeaks', 1);

% Check if any of the dimensions does not have a peak (it is more likely
% that it will be y-direction, but check both). If there is no peak this
% means that there is not a sufficient bearing change to be able to obtain
% localization- the code should exit and report.
if any([isempty(height_y),isempty(height_x)])
    xrange_new=[];
    yrange_new=[];
    flag_locnotpossible = true;
else %Otherwise (if there are peaks in each), obtain ranges
    
    % Take 50% from the peak to be my bounds
    ind_xmarg=find(x_marg>height_x*0.5);

    % If peak is extremely narrow then sometimes not both sides from the peak
    % are included- ensure there are at least 3 values in there with the peak
    % value in the middle:
    if numel(ind_xmarg)<3
        ind_xpeak=find(x_marg==height_x);
        ind_xmarg = [ind_xpeak-1,ind_xpeak,ind_xpeak+1];
    end

    % Take 50% from the peak to be my bounds
    ind_ymarg=find(y_marg>height_y*0.5);

    % If peak is extremely narrow then sometimes not both sides from the peak
    % are included- ensure there are at least 3 values in there with the peak
    % value in the middle (although very unlikely in the y-direction)
    if numel(ind_ymarg)<3
        ind_ypeak=find(y_marg==height_y);
        ind_ymarg = [ind_ypeak-1,ind_ypeak,ind_ypeak+1];
    end

    xrange_new=round([gridparams.x(ind_xmarg(1))-gridparams.dx,gridparams.x(ind_xmarg(end))+gridparams.dx]);
    yrange_new=round([gridparams.y(ind_ymarg(1))-gridparams.dy,gridparams.y(ind_ymarg(end))+gridparams.dy]);
    flag_locnotpossible = false;
end

end