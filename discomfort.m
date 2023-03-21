function [discomfort] = discomfort(cap, SoCinit, chargEff, dischEff, Sobj, Smax, S)
%DISCOMFORT Quantification of discomfort of an EV based on preferences
%   Detailed explanation goes here
baseline = zeros(3,1);
baseline(1) = SoCinit;
for i=2:len(baseline)
    baseline(i) = baseline(i-1) + chargEff*Emax/cap;
    if baseline(i) == Smax
        baseline(i) = Smax;
    end
end

%DEV of Baseline
dev = cap*sum(baseline-S);

%TARget violation 
tar = cap*(Sobj - S(end));

%MILeage
mil = sum(chargeEff*E + e/dischEff) - cap*(S(end) - initSoC);

discomfort = (dev + tar + mil)/3;
end

