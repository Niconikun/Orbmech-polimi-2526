function a = RepeatingGroundTrack(k, m, omega_E, mu)
% RepeatingGroundTrack is a function that compute the semi major axis a for
% a repeating ground track
% 
% PROTOTYPE
%  a = RepeatingGroundTrack(k, m, omega_E, mu)
%
% INPUT:
%   k [1]            number of satellite revolutions
%   m [1]            number of planet revolutions
%   omega_E [1]      Angular velocity of the planet
%   mu [1]           Gravitational parameter of the planet
%
% OUTPUT:
%   a [1]             Semi major axis
% 
% CONTRIBUTORS:
%   Alice Muscas, Federico Masiero
%
% VERSIONS
% 2025-11-10: First Version
% 2025-11-15: Second Version
%
% -------------------------------------------------------------------------

% 1. Mean motion
n = omega_E * k / m;

% 2. Semi major axis
a = (mu / n^2) ^ (1/3);

end