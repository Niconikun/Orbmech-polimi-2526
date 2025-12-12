function [kep] = car2kep(r, v, mu)
% car2kep Conversion from state vector to keplerian element vector
% 
% PROTOTYPE
%  [kep] = car2kep(r, v, mu)
% 
% INPUT
%   r [3x1]          Vector of position
%   v [3x1]          Vector of velocity
%   mu [1]           Gravitational parameter of the Earth
% 
% OUTPUT
%   kep [6x1]        Keplerian element vector
% 
% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico
% 
% VERSIONS
%   2025-10-22: First version
%   2025-10-14: Second version
%   2025-11-12: Third version
% 
% -------------------------------------------------------------------------

% 1. Check
if length(r) ~= 3
    error("r must be a vector of legth 3")
end
if length(v) ~= 3
    error("v must be a vector of legth 3")
end

% 2. Definition of z_vec, x_vec, toll
z_vec = [0; 0; 1];
x_vec = [1; 0; 0];
toll = 1e-3;

% 3. Normalization of r and v
r_norm = norm(r);
v_norm = norm(v);

% 4. Compute radial velocity v_r
v_r = dot(r, v) / r_norm;

% 5. Angular momentum h
h = cross(r, v);
h_norm = norm(h);

% 6. Inclination i
i = acos(h(3) / h_norm);

% 7. Line of node
N = cross(z_vec, h);
N_norm = norm(N);
if N_norm < toll
    N = x_vec;
    N_norm = 1;
end

% 8. Eccentricity
e_vec = (1/mu) * ( (v_norm^2 - mu / r_norm) .* r - (r_norm*v_r) .* v );
e_norm = norm(e_vec);
e = e_norm;
if e_norm  < toll
    e_vec = x_vec;
    e_norm = 1;
    e = 0;
end

% 9. RAAN
if isequal(N, x_vec)
    Omega = 0;
else
    if N(2) >= 0
        Omega = acos(N(1) / N_norm);
    else
        Omega = 2*pi - acos(N(1) / N_norm);
    end
end

% 10. Argument of the pericentre
if e == 0
    omega = 0;
else
    if e_vec(3) >= 0
          omega = acos( dot(N, e_vec) / (e_norm*N_norm) );
    else
          omega = 2*pi - acos( dot(N, e_vec) / (e_norm*N_norm) );
    end
end

% 11. True anomaly
if v_r > 0 
    theta = acos( dot(e_vec, r) ./ (e_norm*r_norm) );
else
    theta = 2*pi - acos( dot(e_vec, r) ./ (e_norm*r_norm) );
end
if e == 0
    omega = theta;
    theta = 0;
end

% 12. Semi major axis
energy = (v_norm^2 / 2) - (mu / r_norm);
a = - mu / (2 * energy);

% 13. Keplerian element vector
kep = [a; e; i; Omega; omega; theta];

end