function plot_AS(AStotal,AS_params,hyph_pos,boat_pos,est_loc, est_loc_bounds)
%plotAS.m is a function that plots the final ambiguity surface (AStotal),
%boat and array position, with the option to plot localization estimates
%when available
%
% INPUTS:
% - AStotal : final ambiguity surface, an A x B matrix where A is number of
%              y coordinates, and B is number of x coordinates 
% - AS_params : a structure containing parameters for ambiguity surface
%               computation
% - hyph_pos : a NxMxT array, where T is number of time steps, N=number of
%              sensors and M is number of coordinates (e.g. x,y or x,y,z)
% - boat_pos : a N x 2 matrix, where 1st column contains x-coordinates
% (in m) and 2nd column contains y-coordinates (in m) of boat position
% - est_loc : nx2 matrix, where each row indicates x and y coordinates 
%               of the estimated source location (in m). n is number of
%               locations/sources.
% - est_loc_bounds : 1 x 2 cell array, where first cell corresponds to min
%   bounds and second cell to max bounds. Each cell is nx2 matrix, 
%   where each row indicates x and y coordinates of the bounds (in m),
%   and n is number of sources.
%      
%
%
% Pina Gruden, Dec 2022, UH Manoa

if nargin < 5 % no localization or error bounds supplied
est_locX=NaN;
est_locY=NaN;
Xbounds=NaN;
Ybounds=NaN;
n=3; % for legend display
elseif nargin < 6 % no error bounds supplied
Xbounds=NaN;
Ybounds=NaN;
  %localization
est_locX=est_loc(:,1)';
est_locY=est_loc(:,2)';
n=4; % for legend display
else % localization and error bounds supplied
    %localization
est_locX=est_loc(:,1)';
est_locY=est_loc(:,2)';
    % error bounds:
xmin=est_loc_bounds{1}(:,1); xmax=est_loc_bounds{2}(:,1);
Xbounds=[xmin,xmax]';
ymin=est_loc_bounds{1}(:,2); ymax=est_loc_bounds{2}(:,2);
Ybounds=[ymin,ymax]';
n=5; % for legend display
end

X= AS_params.X;
Y= AS_params.Y;
xrange=[min(min(AS_params.X)),max(max(AS_params.X))];
yrange=[min(min(AS_params.Y)),max(max(AS_params.Y))];


%figure; 
s{1}=pcolor(X,Y,AStotal);
hold on;
s{1}.EdgeColor='none';
%clim([0,1])
colorbar
axis equal
hph=1;
s{2}=plot(squeeze(hyph_pos(hph,1,:)),squeeze(hyph_pos(hph,2,:)),'r-','Linewidth', 3);
s{3}=plot(boat_pos(:,1), boat_pos(:,2),'g-','Linewidth', 1.5);
s{4}=plot(est_locX,est_locY,'ro', 'MarkerSize',10,'Linewidth', 2);
s{5}=plot(Xbounds,Ybounds,'Color','red','LineStyle','--', 'LineWidth', 2);
s{5}=s{5}';


xlabel(' x (m)'),ylabel('y (m)')
warning('off','MATLAB:legend:IgnoringExtraEntries')
legend([s{1:n}],{'Ambiguity Surface',['Sensor ', num2str(hph),' position'], ...
    'Boat track', 'Estimated source location', 'Estimated perpendicular error bounds'}, ...
    'Location','northoutside')
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
set(gca,'FontSize',16)
hold off

end