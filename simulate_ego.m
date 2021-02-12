function [sim] = simulate_ego(param)
%SIMULATE_EGO
%   INPUTS-
%   param.theta = 0; % facing toward object
%   param.Z = angular variable;
%   param.P = P;
%   param.kappa = 5;
%   param.rp = [75, 75];
%   param.A = 10;
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unpack position
X = param.P(:,2);
Y = param.P(:,3);
T = param.P(:,1); 
fs = mode(diff(T));

% shift param.theta (so that 0 is heading toward refpt)
% param.theta = mod(param.theta-180,360)-180;

% calculate egocentric bearing
ego = deg2rad(mod(atan2d(param.rp(2)-Y, param.rp(1)-X) - param.Z, 360)-180); 

% angular bins
[~, edges, bin] = histcounts(ego, linspace(-pi,pi,101)); % circular?
ctrs = (diff(edges)/2) + edges(1:end-1);
bin(bin==0) = nan; 

% make the von-Mises distribution
[vm_pdf, ~] = circ_vmpdf(ctrs, param.theta, param.kappa);

% lambda matrix
vm_pdf = (param.A.*vm_pdf).*fs;

% loop through each timepoint
simTrn = zeros(length(T), 1);
simTimes = [];
for t = 1:length(T)
    tnow = T(t);
    if isnan(bin(t))
        simTrn(t) = 0;
    else
        lambda = vm_pdf(bin(t));
         n = poissrnd(lambda);
         if n > 0
            simTrn(t) = n;
            simTimes = [simTimes; repmat(tnow,n,1)];
         end
        
    end
end
sim.TR = simTrn;
sim.ST = simTimes;

end

