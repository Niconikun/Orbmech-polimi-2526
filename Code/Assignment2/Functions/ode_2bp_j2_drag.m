function dy = ode_2bp_j2_drag( ~, y, mu, J2, R, omega, c_d, AoverMass)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% J2 is the harmonic coefficient of the planet or object
% R is the radius of the planet being orbited
% omega is the mean angular velocity of the planet being orbited
% c_d is the drag coefficient of the spacecraft
% AoverMass is the area over mass ratio of the spacecraft
% rho is the atmosphere at the position provided
%
%
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-12: First version
% 
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);

altitude = rnorm - R;

x_position = r(1);
y_position = r(2);
z_position = r(3);

v_rel = v - cross(omega,r);

%density
rho = ExponentialaAtmosphericModel(altitude).*1e9; %kg/km^3

% Calculate the drag force
a_drag = -0.5 .* c_d .* AoverMass .* rho .* norm(v_rel).* v_rel;


factor = (3*J2*mu*(R^2))/(2*(rnorm^4));
a_j2 = factor.*[(x_position/rnorm).*(5*((z_position/rnorm)^2) - 1);
    (y_position/rnorm).*(5.*((z_position./rnorm).^2) - 1);
    (z_position./rnorm).*(5.*((z_position./rnorm).^2) - 3)];



% Set the derivatives of the state
dy = [ v;
    ((-mu/rnorm^3)*r)+a_j2 + a_drag];
end