%Rotate xy coordinates to align with x-axis (insert this in localize_track)

p=polyfit(boat_pos(:,1),boat_pos(:,2),1);

phi=atan(p(1)); %the angle between the line and x-axis is inverse tan of the slope

rotmatfnc= @(x) [cos(x),-sin(x);sin(x),cos(x)]; %rotation matrix

boat_pos_rotd=(rotmatfnc(-phi)*boat_pos')';

hyph_pos_rotd= pagetranspose(pagemtimes(rotmatfnc(-phi),pagetranspose(hyph_pos(:,:,:))));