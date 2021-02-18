function [sim] = simulate_hd(param)
%SIMULATE_EGO
% param.theta = 0; % facing toward object
% param.Z = angular variable (range from 0 to 360)
% param.P = P;
% param.kappa = 5;
% param.rp = [75, 75];
% param.A = 10;

% get position
X = param.P(:,2);
Y = param.P(:,3);
T = param.P(:,1); 
fs = mode(diff(T));
Z = param.Z;
pref_theta = mod(param.theta + 90, 360);

% go from 0-360 to -180 to 180 --> & convert to rad
% https://confluence.ecmwf.int/display/CUSF/Longitude+conversion+0~360+to+-180~180
% Z = deg2rad(mod((Z+180),360)-180);
Z = deg2rad(Z);

% angular bins
[~, edges, bin] = histcounts(Z, linspace(0,2*pi,101)); % circular?
ctrs = (diff(edges)/2) + edges(1:end-1);
bin(bin==0) = nan; 

% make the von-Mises distribution
[vm_pdf, ~] = circ_vmpdf(ctrs, deg2rad(pref_theta), param.kappa);

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