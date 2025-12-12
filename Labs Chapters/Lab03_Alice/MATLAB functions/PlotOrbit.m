function [t, S] = PlotOrbit(S0, t, mu, J2, Line, options)
% Plot orbit is a function that propagate and plot an orbit
% 
% PROTOTYPE
%  a = RepeatingGroundTrack(k, m, omega_E, mu)
%
% INPUT:
%   S0 [6x1]         Initial state vector
%   t [1xN]          Time vector
%   mu [1]           Gravitational parameter of the planet
%   J2 [1]           Gravitatonal field constant of the Earth. Equal to zero for the unperturbed case
%   options [-]      Options for the integration
%   Line [-]         Options for the style of the plot
%
% OUTPUT:
%   t [1xM]          Time vector
%   S [6xM]          States vector of the orbit
% 
% CONTRIBUTORS:
%   Alice Muscas
%
% VERSIONS
% 2025-11-24: First Version
%
% -------------------------------------------------------------------------

% Integration
if nargin < 6
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
end

[t, S] = ode45(@(t, y) ode_2bp_perturbed(t, y, mu, J2), t, S0, options);

% Plot
hold on;
plot3(S(:, 1), S(:, 2), S(:, 3), Line, 'LineWidth', 2);
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
grid on; 
axis equal;
title('Orbit Propagation using Two-Body Problem');

end

function dydt = ode_2bp_perturbed(~, y, mu, J2)

    % Ensure y is a column vector
    y = y(:);

    % Extract position and velocity
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    x_r = r(1); 
    y_r = r(2); 
    z_r = r(3);

    % Earth radius
    Re = astroConstants(23);

    % Gravitational acceleration (2-body)
    a_grav = -(mu / r_norm^3) * r;

    % J2 perturbation
    if J2 ~= 0
        factor = 1.5 * J2 * mu * (Re^2) / (r_norm^5);
        zx_ratio = (5 * z_r^2 / r_norm^2);
        aJ2x = factor * x_r * (zx_ratio - 1);
        aJ2y = factor * y_r * (zx_ratio - 1);
        aJ2z = factor * z_r * (zx_ratio - 3);
        aJ2 = [aJ2x; aJ2y; aJ2z];
    else
        aJ2 = [0; 0; 0];
    end

    % Total acceleration
    a_total = a_grav + aJ2;

    % Return derivative
    dydt = [v; a_total];

end