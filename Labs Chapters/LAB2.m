clc
clear
close all
addpath("Functions\");

%% FUNCTIONS
function dy = ode_2bp( ~, y, mu )
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
    % Set the derivatives of the state
    dy = [ v
        (-mu/rnorm^3)*r ];
end

%% LAB CHAPTER 2
% 2 body orbital mechanics

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
tspan = linspace( 0, 10*T, 5000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%for variable compatibility
t_vec = tspan;

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0);
% Plot the results
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;