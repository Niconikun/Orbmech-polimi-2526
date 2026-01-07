function [kep] = car2kep(r, v, mu)
% car2kep Conversion from state vector to keplerian element vector
% 
% Function that compute keplerial element vector from the state vector
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
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-10-22: First version
%   2025-10-14: Second version
%   2025-11-12: Third version
%   2026-01-03: Fourth version


% 0. Definition of z_vec, x_vec, toll
z_vec = [0; 0; 1];
x_vec = [1; 0; 0];
toll = 1e-3;

% 1. Normalization of r and v
r_norm = norm(r);
v_norm = norm(v);

% 2. Compute radial velocity v_r
v_r = dot(r, v) / r_norm;

% 3. Angular momentum h
h = cross(r, v);
h_norm = norm(h);

% 4. Inclination i
i = acos(h(3) / h_norm);

% 5. Line of node
N = cross(z_vec, h);
N_norm = norm(N);
if N_norm < toll
    N = x_vec;
    N_norm = 1;
end

% 6. Eccentricity
e_vec = (1/mu) * ( (v_norm^2 - mu / r_norm) .* r - (r_norm*v_r) .* v );
e_norm = norm(e_vec);
e = e_norm;
if e_norm  < toll
    e_vec = x_vec;
    e_norm = 1;
    e = 0;
end

% 7. RAAN
if isequal(N, x_vec)
    Omega = 0;
else
    if N(2) >= 0
        Omega = real(acos(N(1) / N_norm));
    else
        Omega = real(2*pi - acos(N(1) / N_norm));
    end
end

% 8. Argument of the pericentre
if (i ~= 0) && (i ~= pi)
    if e == 0
        omega = 0;
    else
        if e_vec(3) >= 0
            omega = real(acos( dot(N, e_vec) / (e_norm*N_norm) ));
        else
            omega = real(2*pi - acos( dot(N, e_vec) / (e_norm*N_norm) ));
        end
    end

elseif i == 0
    if e_vec(2) >= 0
          omega = real(acos( dot(N, e_vec) / (e_norm*N_norm) ));
    else
          omega = real(2*pi - acos( dot(N, e_vec) / (e_norm*N_norm) ));
    end

elseif i == pi
    if e_vec(2) <= 0
          omega = real(acos( dot(N, e_vec) / (e_norm*N_norm) ));
    else
          omega = real(2*pi - acos( dot(N, e_vec) / (e_norm*N_norm) ));
    end
end

% 8. True anomaly
if v_r >= -1e-14 
    theta = real(acos( dot(e_vec, r) ./ (e_norm*r_norm) ));
else
    theta = real(2*pi - acos( dot(e_vec, r) ./ (e_norm*r_norm) ));
end
if e == 0
    omega = theta;
    theta = 0;
end

% 9. Semi major axis
energy = (v_norm^2 / 2) - (mu / r_norm);
a = - mu / (2 * energy);

% 10. kep
kep = [a; e; i; Omega; omega; theta];

end