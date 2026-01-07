function [a, T_nodal, Omega_dot] = computeRGTSemiMajorAxis(mu, R, J2, omega, N, M, i, e, tol)
% COMPUTERGTSEMIMAJORAXIS Calculate semi-major axis for J2-perturbed repeating ground track
%
% PROTOTYPE
%   [a, T_nodal, Omega_dot] = computeRGTSemiMajorAxis(mu, R, J2, omega, N, M, i, e, tol)
%
% Inputs:
%   N     - Integer number of orbital revolutions per repeat cycle
%   M     - Integer number of Earth rotations per repeat cycle
%   i     - Inclination [rad]
%   e     - Eccentricity (optional, default = 0)
%   tol   - Tolerance for convergence (optional, default = 1e-9)
%
% Outputs:
%   a         - Required semi-major axis [m]
%   T_nodal   - Nodal (draconic) period [s]
%   Omega_dot - Nodal regression rate [rad/s]
%
% Example:
%   % Calculate for 15 revs per day, 53° inclination
%   [a, T_nodal, Omega_dot] = computeRGTSemiMajorAxis(15, 1, deg2rad(53));
%   altitude = (a - 6378137)/1000; % Altitude in km
%
% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-12: First version
% 
% -------------------------------------------------------------------------

% Physical constants (WGS84)

T_earth = 2*pi/omega;        % Sidereal day [s]

% Default values
if nargin < 8
    e = 0;  % Assume circular orbit
end
if nargin < 9
    tol = 1e-9;  % Convergence tolerance
end

% 1. Keplerian initial guess (simplified repeat condition)
a_kep = (mu * (N * T_earth / (2 * pi * M))^2)^(1/3);

% 2. Iterative solution for J2-compensated semi-major axis
a = a_kep;  % Start with Keplerian guess
max_iter = 100;

for iter = 1:max_iter
    % Current mean motion
    n = sqrt(mu / a^3);
    
    % J2 secular rates (Vallado 4th Ed, Eq. (9-33) to (9-35))
    % Common factor: K = (3/2) * J2 * R^2 * sqrt(mu) / (a^(7/2))
    K = (3/2) * J2 * R^2 * sqrt(mu) / (a^(7/2));
    
    % Nodal regression rate (negative for prograde orbits)
    Omega_dot = -K * cos(i) / (1 - e^2)^2;
    
    % Apsidal rotation rate (not directly needed but included for completeness)
    %omega_dot = K * (2 - (5/2)*sin(i)^2) / (1 - e^2)^2;
    
    % Nodal (draconic) period including J2 effects
    % T_nodal = 2π / (n + ω_dot) for the anomalistic period, but we need:
    % For ground track, we use the nodal period which accounts for Ω_dot
    % The mean motion relative to the ascending node is: n_rel = n - Ω_dot * cos(i)?
    % Actually, the condition is simpler using the repeat formula:
    
    % Calculate the repeat condition error
    T_nodal = 2*pi / n * (1 - (3/2)*J2*(R/a)^2*(1 - (3/2)*sin(i)^2)/sqrt(1-e^2)^3);
    
    % Repeat condition: (ω_earth - Ω̇) * N * T_nodal = 2π * M
    LHS = (omega - Omega_dot) * N * T_nodal;
    RHS = 2 * pi * M;
    
    error = LHS - RHS;
    
    % Check convergence
    if abs(error) < tol
        break;
    end
    
    % Update semi-major axis using Newton-Raphson
    % We need the derivative of error with respect to a
    % Let's use a finite difference approximation
    da = a * 1e-6;  % Small perturbation
    a_temp = a + da;
    
    % Recalculate for perturbed a
    n_temp = sqrt(mu / a_temp^3);
    K_temp = (3/2) * J2 * R^2 * sqrt(mu) / (a_temp^(7/2));
    Omega_dot_temp = -K_temp * cos(i) / (1 - e^2)^2;
    T_nodal_temp = 2*pi / n_temp * (1 - (3/2)*J2*(R/a_temp)^2*(1 - (3/2)*sin(i)^2)/sqrt(1-e^2)^3);
    LHS_temp = (omega - Omega_dot_temp) * N * T_nodal_temp;
    error_temp = LHS_temp - RHS;
    
    % Derivative
    d_error_da = (error_temp - error) / da;
    
    % Newton-Raphson update (with damping factor for stability)
    a = a - 0.5 * error / d_error_da;
    
    % Safety check: keep a within reasonable bounds
    if a < R + 150000  % Below 150 km altitude
        a = R + 150000;
    elseif a > R + 1000000  % Above 10000 km altitude
        a = R + 1000000;
    end
    
    if iter == max_iter
        warning('Maximum iterations reached. Solution may not be fully converged.');
    end
end

% Final nodal period
T_nodal = 2*pi / sqrt(mu/a^3) * (1 - (3/2)*J2*(R/a)^2*(1 - (3/2)*sin(i)^2)/sqrt(1-e^2)^3);
end