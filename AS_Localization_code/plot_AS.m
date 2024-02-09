function plot_AS(AStotal,Gridparams,hyph_pos,boat_pos,est_loc, est_loc_bounds)
%plotAS.m is a function that plots the final ambiguity surface (AStotal),
%boat and array position, with the option to plot localization estimates
%when available
%
% INPUTS:
% - AStotal :  ambiguity surface for each group. It is N x 1 cell where N
% is number of groups. Each cell is an A x B matrix where A is number of
%              y coordinates, and B is number of x coordinates 
% - Gridparams : Grid parameters (X and Y) of each ambiguity surface. It is 
%              N x 1 cell where N is number of groups. Each cell is a 
%              structure containing a matrix of X and Y coordinates.
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

ngroups=size(AStotal,1);

xrange_temp=zeros(ngroups,2);
yrange_temp=zeros(ngroups,2);

% Plot Ambiguity surfaces for each group
for k=1:ngroups

    X= Gridparams{k}.X;
    Y= Gridparams{k}.Y;
    xrange_temp(k,:)=[min(min(Gridparams{k}.X)),max(max(Gridparams{k}.X))];
    yrange_temp(k,:)=[min(min(Gridparams{k}.Y)),max(max(Gridparams{k}.Y))];

    s{1}=pcolor(X,Y,AStotal{k});
    hold on;
    s{1}.EdgeColor='none';
    colorbar
    axis equal

end

% plot boat track, hydrophone positions, estimated locations and bounds
hph=1;
hyph_pos_x=squeeze(hyph_pos(hph,1,:));
hyph_pos_y=squeeze(hyph_pos(hph,2,:));
s{2}=plot(hyph_pos_x,hyph_pos_y,'r-','Linewidth', 3);
s{3}=plot(boat_pos(:,1), boat_pos(:,2),'g-','Linewidth', 1.5);
s{4}=plot(est_locX,est_locY,'ro', 'MarkerSize',10,'Linewidth', 2);
s{5}=plot(Xbounds,Ybounds,'Color','red','LineStyle','--', 'LineWidth', 2);
s{5}=s{5}';

%compute x and y coordinate range limits:
xrange_AS=[min(xrange_temp(:,1)),max(xrange_temp(:,2))]; %min and max bounds in x-coordinate for AS
yrange_AS=[min(yrange_temp(:,1)),max(yrange_temp(:,2))]; %min and max bounds in y-coordinate for AS

xrange = [min([xrange_AS(1),min(boat_pos(:,1)),min(hyph_pos_x)]), ...
    max([xrange_AS(2),max(boat_pos(:,1)),max(hyph_pos_x)])];
yrange = [min([yrange_AS(1),min(boat_pos(:,2)),min(hyph_pos_y)]), ...
    max([yrange_AS(2),max(boat_pos(:,2)),max(hyph_pos_y)])];

xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
xlabel(' x (m)'),ylabel('y (m)')
warning('off','MATLAB:legend:IgnoringExtraEntries')
legend([s{1:n}],{'Ambiguity Surface',['Sensor ', num2str(hph),' position'], ...
    'Boat track', 'Estimated source location', 'Estimated perpendicular error bounds'}, ...
    'Location','northoutside')
set(gca,'FontSize',16)
hold off

end