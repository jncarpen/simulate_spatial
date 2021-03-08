function [sim] = simulate_placehd(map, param)
%SIMULATE_PLACEMOD
%   param.P:        position vector
%   param.theta:    head direction (range 0-360 deg)
%   param.Z:        head direction for each timestamp (0 should be 'up')
%   param.kappa:    steepness of von-mises distribution (smaller k = wider);
%                   a large value (i.e. 10) will give a purely directional
%                   place cell. 
%   param.A:        peak firing rate of unit
%   Jordan Carpenter. 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort inputs
t = param.P(:,1);
Z = deg2rad(param.Z);
pref_theta = mod(param.theta + 90, 360);
pdx = map.pdx;
tpf = mode(diff(t)); % time per frame

% angular bins
[~, edges, bin] = histcounts(Z, linspace(0,2*pi,101)); % circular?
ctrs = (diff(edges)/2) + edges(1:end-1);
bin(bin==0) = nan; 

% make the von-Mises distribution
[vm_pdf, ~] = circ_vmpdf(ctrs, deg2rad(pref_theta), param.kappa);

simTrn = zeros(length(t), 1);
simTimes = [];
for frame = 1:length(t)
    
    % timestamp now
    tnow = t(frame);
    
    % find x,y bin for current sample
    xnow = map.whichBin.x(frame);
    ynow = map.whichBin.y(frame);
    znow = bin(frame);
    
    % if timestamp has a nan
    if isnan(xnow) | isnan(ynow) | isnan(znow)
        simTrn(frame) = 0;
    else
        
        % find probability (von-mises)- hd
        probZ = vm_pdf(znow);
        
        % find probability (gaussian)- spatial
        probSpatial = pdx(xnow, ynow);
        
        % multiply the probabilities together & scale by peak fr
        lambda = (param.A.*(probZ.*probSpatial)).*tpf;
    
        %draw random number of spikes (n) from Poisson distribution
        n = poissrnd(lambda);
    
        %store number of spikes and spike-to-position indices
        if n>0
            simTrn(frame) = n;
            simTimes = [simTimes; repmat(tnow,n,1)];
        end 
    
    end
end

sim.TR = simTrn;
sim.ST = simTimes;

end
