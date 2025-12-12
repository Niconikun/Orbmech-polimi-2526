function [dv_tot, dv1, dv2] = VelocityCost(r1, v1, r2, v2, ToF, mu)
% VelocityCost is a function that compute total velocity cost for a mission
% 
% PROTOTYPE
%  [dv_tot, dv1, dv2] = VelocityCost(r1, v1, r2, v2, ToF, mu)
%
% INPUT:
%   r1 [1x3]         Initial position
%   v1 [1x3]         Initial velocity
%   r2 [1x3]         Final position
%   v2 [1x3]         Final velocity
%   ToF [1]          Time of flight in seconds
%   mu [1]           Planetary constant
%
% OUTPUT:
%   dv_tot [1]       Total velocity cost for the transfer maneuver
%   dv1 [1]          Velocity cost for the first maneuver
%   dv2 [1]          Velocity cost for the second maneuver
% 
% CONTRIBUTORS:
%   Alice Muscas
%
% VERSIONS
% 2025-11-23: First Version
%
% -------------------------------------------------------------------------

[~, ~, ~, ~, v1_T, v2_T, ~, ~] = lambertMR(r1, r2, ToF, mu, 0, 0, 0, 0);
v1_T = v1_T';
v2_T = v2_T';

dv1 = norm(v1_T - v1);
dv2 = norm(v2 - v2_T);

dv_tot = dv1 + dv2;

end