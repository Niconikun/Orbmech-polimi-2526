function [r, v] = kep2car(kep, mu)
% kep2car compute the state vector in orbital frame from the keplerian
% element vector
% 
% PROTOTYPE
%  [r, v] = kep2car(kep, mu)
% 
% INPUT
%   kep [6x1]        Keplerian element vector
%   mu [1]           Gravitational parameter of the Earth
% 
% OUTPUT
%   r [3x1]          Position vector
%   v [3x1]          Velocity vector
% 
% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-10-22: First version
%   2025-10-14: Second version
% 
% -------------------------------------------------------------------------

% 1. Extract component from keplerian element vector
a = kep(1);
e = kep(2);
i = kep(3);
Omega = kep(4);
omega = kep(5);
theta = kep(6);

% 2. Check
if length(kep) ~= 6
    error("Keplerian element vector is not complete")
end

% 3. Angular momentum
p = a * (1 - e^2);
h = sqrt(p * mu);

% 4. r and v in perifocal frame
r = [cos(theta); sin(theta); 0];
r = (h^2 / mu) / (1 + e*cos(theta)) .* r;
v = [-sin(theta); e + cos(theta); 0];
v = (mu / h) .* v;

% 5. Rotation matrix
R1 = [cos(omega)   sin(omega)   0;
      -sin(omega)  cos(omega)   0;
               0                    0           1];
R2 = [ 1     0        0;
       0   cos(i)   sin(i);
       0  -sin(i)   cos(i)];
R3 = [ cos(Omega)  sin(Omega)  0;
      -sin(Omega)  cos(Omega)  0;
           0           0       1];
A = (R1 * R2 * R3)';

% 6. r and v in orbit frame
r = A * r;
v = A * v;

end