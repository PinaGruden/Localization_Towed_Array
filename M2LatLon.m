function PosLatLon = M2LatLon(PosM, PosZero)

% Convert position in Meters to Lat/Lon according to a given zero position (in
% lat/lon)
% PosLatLon = [lat(decimal degrees)  lon(decimal degrees)];
% PosZero = [lat(decimal degrees)  lon(decimal degrees)];
% PosMeters = [x y] in m

dlat = conv_m2lat(PosM(:,2), PosZero(1));
dlon = conv_m2lon(PosM(:,1), PosZero(1));

lat = PosZero(1) + dlat;
lon = PosZero(2) + dlon;

PosLatLon = [lat lon];

if size(PosM,2) == 3;
    PosLatLon = [PosLatLon PosM(:,3)];
end

function dlon = conv_m2lon(dx, alat)
%
% dlon = longitude difference in degrees
% dx   = longitude difference in meters
% alat = average latitude between the two fixes

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dlon = dx ./ p;

function dlat = conv_m2lat(dy, alat)
%
% dy   = latitude difference in meters
% dlat = latidute difference in degrees
% alat = average latitude between the two fixes
% Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 

rlat = alat * pi/180;
m = 111132.09 * ones(size(rlat)) - ...
    566.05 * cos(2 * rlat) + 1.2 * cos(4 * rlat);
dlat = dy ./ m;