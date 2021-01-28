function [sim] = simulate_place(map, t)
%SIMULATE_PLACE_V2 
%   INPUTS -
%   t:                    Tx1 vector of timestamps (s)
%   map.z:                ratemap (Hz)
%   map.whichBin.x:       bin the animal occupies for each frame (x)
%   map.whichBin.y:       bin the animal occupies for each frame (y)
%   OUTPUTS -
%   

%   J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ratemap
RM = map.z;

% sampling freq (s)
fs = mode(diff(t));

% multiply by sample time to get spikes/frame (or spikes/position sample)
lambdaMatrix = RM.*fs;

simTrn = zeros(length(t), 1);
simTimes = [];
for frame = 1:length(t)
    
    % timestamp now
    tnow = t(frame);
    
    % find x,y bin for current sample
    xnow = map.whichBin.x(frame);
    ynow = map.whichBin.y(frame);
    
    % if timestamp has a nan
    if isnan(xnow) || isnan(ynow)
        simTrn(frame) = 0;
        
    else
        
    % find value of ratemap (fr now)
    lambda = lambdaMatrix(xnow, ynow);
    
    %draw random number of spikes (n) from Poisson distribution
    n = poissrnd(lambda);
    
    %store number of spikes and spike-to-position indices
    if n>0
        simTrn(frame) = n;
        simTimes = [simTimes; tnow];
    end 
    
    end
end

sim.TR = simTrn;
sim.ST = simTimes;

end

