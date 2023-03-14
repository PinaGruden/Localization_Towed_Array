function plot_AS(AStotal,AS_params,hyph_pos,boat_pos,est_loc_m)
%plotAS.m is a function that plots the final ambiguity surface (AStotal),
%boat and array position.
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
% - est_loc_m : nx2 matrix, where each row indicates x and y coordinates 
%               of the estimated source location (in m). n is number of
%               locations.
%
%
% Pina Gruden, Dec 2022, UH Manoa

X= AS_params.X;
Y= AS_params.Y;
xrange=[min(min(AS_params.X)),max(max(AS_params.X))];
yrange=[min(min(AS_params.Y)),max(max(AS_params.Y))];
nlocs=size(est_loc_m,1);

%figure; 
s=pcolor(X,Y,AStotal);
hold on;
s.EdgeColor='none';
%clim([0,1])
colorbar
axis equal
hph=1;
plot(squeeze(hyph_pos(hph,1,:)),squeeze(hyph_pos(hph,2,:)),'r-','Linewidth', 3)
plot(boat_pos(:,1), boat_pos(:,2),'g-','Linewidth', 1.5)
for n=1:nlocs
plot(est_loc_m(n,1),est_loc_m(n,2),'ro', 'MarkerSize',10,'Linewidth', 2) 
end

xlabel(' x (m)'),ylabel('y (m)')
legend({'Ambiguity Surface',['Sensor ', num2str(hph),' position'],'Boat track', 'Estimated source location'})
xlim([xrange(1),xrange(2)])
ylim([yrange(1),yrange(2)])
set(gca,'FontSize',16)
hold off

end