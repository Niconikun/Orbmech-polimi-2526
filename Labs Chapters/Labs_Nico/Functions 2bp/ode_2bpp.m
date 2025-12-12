function dy = ode_2bpp( ~, y, mu, J2, R )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);

x_position = r(1);
y_position = r(2);
z_position = r(3);

factor = (3*J2*mu*(R^2))/(2*(rnorm^4));
a_j2 = factor.*[(x_position/rnorm).*(5*((z_position/rnorm)^2) - 1);
    (y_position/rnorm).*(5.*((z_position./rnorm).^2) - 1);
    (z_position./rnorm).*(5.*((z_position./rnorm).^2) - 3)];
% Set the derivatives of the state
dy = [ v;
    ((-mu/rnorm^3)*r)+a_j2];
end