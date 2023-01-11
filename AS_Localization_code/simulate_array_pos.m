function [hyph_pos,x_boat, y_boat] = simulate_array_pos(Ntsteps,params)
% simulate_array_pos.m simulates sensor positions assuming the boat moves
% in a straight line specified by the slope and interstect
%
% INPUTS:
% - Ntsteps : number of time steps to simulate
% - params : a structure that contains information on the boat and line
%            along which the boat moves- the fileds are:
%           ~ params.mb = line gradient
%           ~ params.b = intersect with y axis
%           ~ params.boatspeed = boat speed in m/s
%           ~ params.timestep = how much time elapses in each time step in s
%
% OUTPUTS:
% - hyph_pos : a NxMxT array, where T= number of time steps, N=number of
%              sensors and M is number of coordinates (e.g. x,y or x,y,z)
%
%
% Pina Gruden, Dec 2022, UH Manoa

% Simulate Boat movement and Hydrophone positions
        % get boat positions along the specified line:
        x_line=0:1:Ntsteps+1;
        y_line=params.mb.*x_line+params.b;
        dxy_line=sqrt((x_line(2)-x_line(1))^2+(y_line(2)-y_line(1))^2);

        boatmove_m = params.boatspeed*params.timestep; %distance the boat traveles in one time step
        x_boat= x_line(1:end-1) + boatmove_m.*(x_line(2:end)-x_line(1:end-1))./dxy_line;
        y_boat=params.mb.*x_boat+params.b;

        hyph_pos=nan(2,2,Ntsteps);%[x1,y1; x2,y2];
        %first sensor is at the position of the boat
        hyph_pos(1,1,:)=x_boat(1:end-1);
        hyph_pos(1,2,:)=y_boat(1:end-1);
        %second sensor is distance d behind first
        %hyph_pos(2,1,:)=hyph_pos(1,1,1:end-1) - d.*(hyph_pos(1,1,2:end)-hyph_pos(1,1,1:end-1))./dxy_line;
        hyph_pos(2,1,:)=x_boat(1:end-1) - params.d.*(x_boat(2:end)-x_boat(1:end-1))./dxy_line;
        hyph_pos(2,2,:)=params.mb.*hyph_pos(2,1,:)+params.b;


end