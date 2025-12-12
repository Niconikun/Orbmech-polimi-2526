function epsilon = specificenergy(v, mu, r, a)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    if ~isempty(v) && ~isempty(r) && ~isempty(mu)
        epsilon = (v.^2)/2 - mu/r; % Calculate specific energy
    elseif ~isempty(mu) && ~isempty(a)
        epsilon = -mu/(2.*a); % Assign epsilon to a if mu and a are provided
    else
        epsilon = NaN; % Return NaN if insufficient inputs are provided
    end
end