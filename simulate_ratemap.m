function [map] = simulate_ratemap(root)
%SIMULATE_RATEMAP Summary of this function goes here
%   INPUTS -
%   root.peakrate:      peak rate (in Hz)
%   root.ctrX:          center of place field (x)          
%   root.ctrY:          center of place field (y)
%   root.sigmaX:        spread of place field (x)
%   root.sigmaY:        spread of place field (y)
%   root.size:          size of square arena
%   root.bins:          # of bins (in each dimension)
%   root.P:             position vector, optional
%   OUTPUTS - 
%   map.z:              ratemap (Hz)
%   map.whichBin:       bin number for each timestamp
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude of gaussian (max FR in Hz)
A = root.peakrate;

% position of gaussian center
x0 = root.ctr(1);
y0 = root.ctr(2);

% size of arena
sz = root.size;
bins = root.bins;
bin_sz = sz/bins; % in cm

% sigma
sigmaX = root.sigma(1);
sigmaY = root.sigma(2);

% bin arena (square)
xEdges = linspace(0,sz,bins+1);
yEdges = linspace(0,sz,bins+1);
xCenter = (diff(xEdges)./2)+xEdges(1:end-1);
yCenter = (diff(yEdges)./2)+yEdges(1:end-1);
[X,Y] = meshgrid(xCenter, yCenter);

xpos = root.P(:,2); ypos = root.P(:,3);
[~, ~, xbin] = histcounts(xpos, xEdges);
xbin(xbin == 0) = nan;  % x=0 indicates values outside edges range

[~, ~, ybin] = histcounts(ypos, yEdges);
ybin(ybin == 0) = nan;  % x=0 indicates values outside edges range 

% generate ratemap (Hz)
fxy = A.*exp(-(((X-x0).^2)./(2*sigmaX.^2) + ((Y-y0).^2.)/(2*sigmaY.^2)));

% output
map.z = fxy;
map.whichBin.x = xbin;
map.whichBin.y = ybin;

end
