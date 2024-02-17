function [gridparams,wpos2D] = make_2Dgrids(xrange,yrange,dx,dy)
%function make_2Dgrids.m creates a 2-D grid and also outputs a structure of
%grid parameters
%
% INPUTS:
% -xrange - 1 x 2 vector of range [min, max] for x-coordinate
% -yrange - 1 x 2 vector of range [min, max] for y-coordinate
% -dx - resolution in m in x direction
% -dy - resolution in m in y direction
%
% OUTPUTS:
% - gridparams : a structure contining grid parameters, such as x&y-ranges,
%               resolution, number of grid points
% - wpos2D : a N x 2 matrix, where N is number of all grid points and the
%           columns correspond to x and y coordinates respectively
%
% Pina Gruden, Jan 2024, UH Manoa


% check that x & y range are not just one number, otherwise add 200 m to
% the range (100 m to each side)
if diff(xrange)==0
xrange= [xrange(1)-100; xrange(2)+100];
end

if diff(yrange)==0
yrange= [yrange(1)-100; yrange(2)+100];
end

% - Grid range and resolution:

gridparams.xrange=xrange; % x range in m
gridparams.yrange=yrange; % y range in m
gridparams.dx=dx; % grid step size in m (resolution) in x direction
gridparams.dy=dy;% grid step size in m (resolution) in y direction
gridparams.x= xrange(1):gridparams.dx:xrange(2);
gridparams.Ngp_x=length(gridparams.x);
gridparams.y= yrange(1):gridparams.dy:yrange(2);
gridparams.Ngp_y=length(gridparams.y);

[X,Y] = meshgrid(gridparams.x,gridparams.y);
gridparams.X=X(:,:,1);
gridparams.Y=Y(:,:,1);
wpos2D=[X(:),Y(:)]; %gives grid of N x [x,y] coordinates

end