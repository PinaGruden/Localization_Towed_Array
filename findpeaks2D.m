function [peakdata] =findpeaks2D(X,Y,LStotal)

%Find peaks in the ambiguity surface and compute their widths

% Find dimensions to set up loop
xdim = size(X,1);
ydim = size(Y,2);

% Loop through x dimension to find peaks of each row
xpeaks = zeros(size(LStotal));
xwidths = NaN(size(LStotal));
for i = 1:xdim
    [~,locs,w] = findpeaks(LStotal(i,:));
    xpeaks(i,locs) = 1;
    xwidths(i,locs) = w;
end

% Loop through y dimension to find peaks of each row
ypeaks = zeros(size(LStotal));
ywidths = NaN(size(LStotal));
for i = 1:ydim
    [~,locs,w] = findpeaks(LStotal(:,i));
    ypeaks(locs,i) = 1;
    ywidths(locs,i) = w;
end

% Find indices that were peaks in both x and y
peak_inds = xpeaks+ypeaks == 2;

% Plot
h=figure;
hold on
% peaks
s=surf(X,Y,LStotal);
s.EdgeColor='none';
colorbar
%axis equal
plot3(X(peak_inds),Y(peak_inds),LStotal(peak_inds),'r*','MarkerSize',12,'Linewidth',2)
xlabel('x (m)'),ylabel('y (m)'), zlabel('Ambiguity value')
fontsize(h,14,'points')

% Save data to sruct
peakdata = struct;
peakdata.peakZ = LStotal(peak_inds);
peakdata.peakX = X(peak_inds);
peakdata.peakY = Y(peak_inds);
peakdata.peakXWidth = xwidths(peak_inds);
peakdata.peakYWidth = ywidths(peak_inds);

end