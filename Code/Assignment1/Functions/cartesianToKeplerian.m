function [a, e, i, omega, RAAN, theta] = cartesianToKeplerian(x, y, z, vx, vy, vz)
% Converts Cartesian coordinates to Keplerian parameters
% Inputs:
%   x, y, z - Position coordinates (in kilometers)
%   vx, vy, vz - Velocity coordinates (in kilometers per second)
% Outputs:
%   a - Semi-major axis (in kilometers)
%   e - Eccentricity
%   i - Inclination (in degrees)
%   omega - Argument of periapsis (in degrees)
%   RAAN - Right Ascension of Ascending Node (in degrees)
%   theta - True anomaly (in degrees)

% Gravitational parameter for Earth (km^3/s^2)
mu = 398600;

% Position vector
r = [x; y; z];
% Velocity vector
v = [vx; vy; vz];

% Specific angular momentum
h = cross(r, v);
h_mag = norm(h);

% Eccentricity vector
e_vec = (cross(v, h) / mu) - (r / norm(r));
e = norm(e_vec);

% Semi-major axis
energy = norm(v)^2 / 2 - mu / norm(r);
a = -mu / (2 * energy);

% Inclination
i = acos(h(3) / h_mag) * (180 / pi);

% Right Ascension of Ascending Node
RAAN = atan2(h(1), -h(2)) * (180 / pi);
if RAAN < 0
    RAAN = RAAN + 360;
end

% Argument of periapsis
nu = atan2(e_vec(3), dot(e_vec, r) / norm(r));
omega = atan2(e_vec(2), e_vec(1)) * (180 / pi);
if omega < 0
    omega = omega + 360;
end

% True anomaly
theta = atan2(dot(r, cross(e_vec, h)) / (norm(r) * h_mag), dot(e_vec, r) / norm(r)) * (180 / pi);
if theta < 0
    theta = theta + 360;
end
