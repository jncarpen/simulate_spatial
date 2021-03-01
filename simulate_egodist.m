function [sim] = simulate_egodist(param)
%SIMULATE_EGODIST
%   INPUTS-
%   param.A = firing rate (?); default=100 to 300;
%   param.theta = 0; % facing toward object
%   param.Z = angular variable;
%   param.P = P;
%   param.rp = [75, 75];
%   param.sigma = 8; --> larger sigma = wider rings
%   param.radius = 25;
%   param.kappa = 3;
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort inputs
X = param.P(:,2);
Y = param.P(:,3);
t = param.P(:,1);
rp = param.rp;
sigma = param.sigma;
Z = deg2rad(param.Z);
tpf = mode(diff(t)); % time per frame

% shift theta
pref_theta = mod(param.theta + 90, 360);
pref_theta = deg2rad(pref_theta);

% calculate egocentric bearing
allo = mod(atan2d(param.rp(2)-Y, param.rp(1)-X), 360);
allo = mod(allo + 90, 360);
ego = deg2rad(mod(allo - param.Z, 360));

% angular bins
[~, edges, bin] = histcounts(ego, linspace(0,2*pi,101)); % circular?
ctrs = (diff(edges)/2) + edges(1:end-1);
bin(bin==0) = nan; 

% make the von-Mises distribution
[vm_pdf, ~] = circ_vmpdf(ctrs, pref_theta, param.kappa);

% gaussian distribution for distance
d= sqrt((rp(1) - X).^2+ (rp(2) - Y).^2);
[~, dedges, dbin] = histcounts(d, linspace(0,nanmax(d),101));
dctrs = (diff(dedges)/2) + dedges(1:end-1);
dfxy = (1/(sigma.*sqrt(2*pi))) .* exp((-(dctrs-param.radius).^2)./(2.*(sigma.^2)));

simTrn = zeros(length(t), 1);
simTimes = [];
for frame = 1:length(t)
    
    % timestamp now
    tnow = t(frame);
    
    % find x,y bin for current sample
    dnow = dbin(frame);
    znow = bin(frame);
    
    % if timestamp has a nan
    if isnan(dnow) | isnan(znow)
        simTrn(frame) = 0;
    else
        
        % find probability (von-mises)- hd
        probZ = vm_pdf(znow);
        
        % find probability (gaussian)- spatial
        probD = dfxy(dnow);
        
        % multiply the probabilities together & scale by peak fr
        lambda = (param.A.*(probZ.*probD)).*tpf;
    
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

