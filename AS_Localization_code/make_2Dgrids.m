function [gridparams,wpos2D] = make_2Dgrids(xrange,yrange,dx,dy)

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