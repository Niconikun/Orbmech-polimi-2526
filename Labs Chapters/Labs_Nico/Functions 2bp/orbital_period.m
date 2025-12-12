function period = orbital_period(mu, a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    period = 2.*pi.*sqrt((a.^3)./mu);
end