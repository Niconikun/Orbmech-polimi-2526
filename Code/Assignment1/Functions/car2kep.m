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

% 1. Check and ensure vectors are column vectors
if length(r) ~= 3
    error("r must be a vector of length 3")
end
if length(v) ~= 3
    error("v must be a vector of length 3")
end

% Convert to column vectors if necessary
r = r(:);
v = v(:);

% 2. Definition of z_vec, x_vec, toll
z_vec = [0; 0; 1];
x_vec = [1; 0; 0];
toll = 1e-10;  % Reduced tolerance

% 3. Normalization of r and v
r_norm = norm(r);
v_norm = norm(v);

% 4. Compute radial velocity v_r
v_r = dot(r, v) / r_norm;

% 5. Angular momentum h
h = cross(r, v);
h_norm = norm(h);

% Handle special case of zero angular momentum (degenerate orbit)
if h_norm < toll
    % Use a default direction for h
    h = [0; 0; 1];
    h_norm = 1;
end

% 6. Inclination i
i = acos(h(3) / h_norm);

% 7. Line of node
N = cross(z_vec, h);
N_norm = norm(N);
if N_norm < toll
    N = x_vec;
    N_norm = 1;
end

% Ensure N is a column vector
N = N(:);

% 8. Eccentricity vector
e_vec = (1/mu) * ((v_norm^2 - mu/r_norm) * r - (r_norm * v_r) * v);
e_norm = norm(e_vec);
e = e_norm;

if e_norm < toll
    e_vec = x_vec;  % Use default direction for circular orbits
    e_norm = 1;
    e = 0;
else
    e_vec = e_vec(:);  % Ensure column vector
end

% 9. RAAN (Right Ascension of Ascending Node)
if isequal(N, x_vec)
    Omega = 0;
else
    cos_Omega = N(1) / N_norm;
    % Clamp to [-1, 1] to avoid numerical issues
    cos_Omega = max(-1, min(1, cos_Omega));
    if N(2) >= 0
        Omega = acos(cos_Omega);
    else
        Omega = 2*pi - acos(cos_Omega);
    end
end

% 10. Argument of the pericentre
if e < toll
    omega = 0;
else
    cos_omega = dot(N, e_vec) / (e_norm * N_norm);
    % Clamp to [-1, 1]
    cos_omega = max(-1, min(1, cos_omega));
    if e_vec(3) >= 0
        omega = acos(cos_omega);
    else
        omega = 2*pi - acos(cos_omega);
    end
end

% 11. True anomaly
if e < toll
    % For circular orbits, use the angle from the node line
    if N_norm > toll
        cos_theta = dot(N, r) / (r_norm * N_norm);
        cos_theta = max(-1, min(1, cos_theta));
        if r(3) >= 0
            theta = acos(cos_theta);
        else
            theta = 2*pi - acos(cos_theta);
        end
    else
        theta = 0;
    end
else
    cos_theta = dot(e_vec, r) / (e_norm * r_norm);
    cos_theta = max(-1, min(1, cos_theta));
    if v_r >= 0
        theta = acos(cos_theta);
    else
        theta = 2*pi - acos(cos_theta);
    end
end

% 12. Semi major axis
energy = (v_norm^2 / 2) - (mu / r_norm);
if abs(energy) < toll
    a = Inf;  % Parabolic orbit
else
    a = -mu / (2 * energy);
end

% 13. Keplerian element vector
kep = [a; e; i; Omega; omega; theta];

% Normalize angles to [0, 2π)
kep(6) = mod(kep(6), 2*pi);  % True anomaly
kep(5) = mod(kep(5), 2*pi);  % Argument of periapsis
kep(4) = mod(kep(4), 2*pi);  % RAAN

end