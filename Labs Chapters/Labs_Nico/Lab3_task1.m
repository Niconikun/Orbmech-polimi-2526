clc
clearvars
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);  % Earth radius [km]
J2_E = astroConstants(33); % Earth J2

% Earth rotation parameters (sidereal rotation)
w_earth = 360 / 86164; % deg/s (sidereal rotation rate)
t_0 = 0;
theta_g_t_0 = 0; % Initial Greenwich sidereal time [deg]

% Initial condition
kep_parameters = [80795.9, 0.1976, 76.7, 305.3107, 267.2574, 0]; % a [km],e [-],i [deg],Omega [deg],omega [deg],theta [deg]

[r0, v0] = kep2car(kep_parameters, mu_E);
y0 = [r0; v0];  % Combine position and velocity into 6-element vector

% Set time span
a = kep_parameters(1); % Semi-major axis [km]
T_period = 2*pi*sqrt(a^3/mu_E); % Orbital period [s]
tspan = linspace( 0, 5*T_period, 10000 ); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan, y0, options );

% Calculate ground track
r = Y(:,1:3);
gtrack_data = zeros(length(T),4);

% In your main script, before the ground track loop:
fprintf('First position vector: [%f, %f, %f]\n', r(1,1), r(1,2), r(1,3));
fprintf('Time range: %f to %f seconds\n', tspan(1), tspan(end));
fprintf('Earth rotation rate: %f deg/s\n', w_earth);

for i = 1:length(T)
    R = norm(r(i,:));
    %fprintf('Position magnitude: %f km\n', R); % Debug
    
    % Calculate declination (latitude)
    gtrack_data(i,1) = asin(r(i,3) / R); %delta
    delta = gtrack_data(i,1);
    %fprintf('Declination (delta): %f rad, %f deg\n', delta, rad2deg(delta));
    
    % Calculate right ascension
    gtrack_data(i,2) = atan2(r(i,2), r(i,1)); %alpha
    alpha = gtrack_data(i,2);
    %fprintf('Right ascension (alpha): %f rad, %f deg\n', alpha, rad2deg(alpha));
    
    % Calculate Greenwich sidereal time
    w_planet_rad = deg2rad(w_earth);
    theta_g = theta_g_t_0 + w_planet_rad * (T(i) - t_0);
    %fprintf('Greenwich sidereal time (theta_g): %f rad, %f deg\n', theta_g, rad2deg(theta_g));
    
    % Calculate longitude
    lon = rad2deg(alpha - theta_g);
    gtrack_data(i,3) = mod(lon + 180, 360) - 180;
    lon = gtrack_data(i,3);
    
    % Calculate latitude
    gtrack_data(i,4) = rad2deg(delta);
    lat = gtrack_data(i,4);
    
    %fprintf('Final output - Lon: %f deg, Lat: %f deg\n\n', lon, lat);
end

textureImage = imread('EarthTexture.jpg');

% Create the sphere
[XS,YS,ZS] = sphere(50);

scaleFactor = 6371; % Radius of Earth in kilometers; adjust as needed

% Scale the sphere
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;

% Plot the results
figure(1)

hold on
plot3( Y(:,1), Y(:,2), Y(:,3))
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none');


% Finalize the plot
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

% Improved ground track plot that handles wrapping
figure(2)
textureImage = imread('EarthTexture.jpg');
imagesc([-180 180], [-90 90], flipud(textureImage));
set(gca, 'YDir', 'normal');
hold on

% Handle longitude wrapping for continuous lines
lon = gtrack_data(:,3);
lat = gtrack_data(:,4);

% Detect and handle discontinuities
wrap_indices = find(abs(diff(lon)) > 180);
starts = [1; wrap_indices + 1];
ends = [wrap_indices; length(lon)];

% Plot each continuous segment
for i = 1:length(starts)
    idx = starts(i):ends(i);
    plot(lon(idx), lat(idx), 'r', 'LineWidth', 1.5);
end

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track')
xlim([-180 180])
ylim([-90 90])
grid on
hold off

%% ephemeris comparison - Fixed and Optimized
clearvars
clc
close all
addpath("Functions 2bp\")

% Constants
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
J2_E = astroConstants(33); % Earth J2
R_E = astroConstants(23);  % Earth radius [km]

% Earth rotation parameters (sidereal rotation)
w_earth_rad = 2*pi / 86164;     % rad/s (sidereal)
w_earth = rad2deg(w_earth_rad);  % deg/s (for consistency)
theta_g_t0_rad = deg2rad(0);    % initial GST in rad
t_0 = 0; % reference time for propagation

%% Parse ephemeris data
[utc_time, mission_info] = read_ephemeris_file('Ephemeris_Chandra.txt');

% Extract Keplerian elements in correct order: [a, e, i, Ω, ω, θ]
% Note: mission_info(:,10) is semi-major axis, mission_info(:,2) is periapsis distance
% We need to compute eccentricity from periapsis distance and semi-major axis
a_eph = mission_info(:,10);      % Semi-major axis [km]
e_eph = mission_info(:,1);       % Eccentricity (already in column 1)
i_eph = deg2rad(mission_info(:,3));  % Inclination [rad]
OM_eph = deg2rad(mission_info(:,4)); % RAAN [rad]
w_eph = deg2rad(mission_info(:,5));  % Argument of periapsis [rad]
TA_eph = deg2rad(mission_info(:,9)); % True anomaly [rad]

keplerian_elements_eph = [a_eph, e_eph, i_eph, OM_eph, w_eph, TA_eph];

%% Process ephemeris ground track
fprintf('Processing ephemeris ground track...\n');
gtrack_data_eph = zeros(length(utc_time), 4);

for i = 1:length(utc_time)
    % Convert to Cartesian coordinates
    [r_eph, v_eph] = kep2car(keplerian_elements_eph(i,:), mu_E);
    
    % Time difference from reference
    dt = seconds(utc_time(i) - utc_time(1));
    R_eph = norm(r_eph);
    
    % Calculate declination (latitude in radians)
    delta_eph = asin(r_eph(3) / R_eph);
    gtrack_data_eph(i,1) = delta_eph;
    
    % Calculate right ascension
    alpha_eph = atan2(r_eph(2), r_eph(1));
    gtrack_data_eph(i,2) = alpha_eph;
    
    % Calculate Greenwich sidereal time
    theta_g_eph = theta_g_t0_rad + w_earth_rad * dt;
    
    % Calculate longitude (wrap to [-180, 180])
    lon_eph = rad2deg(alpha_eph - theta_g_eph);
    gtrack_data_eph(i,3) = mod(lon_eph + 180, 360) - 180;
    
    % Calculate latitude in degrees
    gtrack_data_eph(i,4) = rad2deg(delta_eph);
end

%% Propagate orbit from initial conditions
fprintf('Propagating orbit...\n');

% Initial conditions from first ephemeris point
[r0, v0] = kep2car(keplerian_elements_eph(1,:), mu_E);
y0 = [r0; v0];

% Set time span based on orbital period
a = keplerian_elements_eph(1,1); % Semi-major axis [km]
T_period = 2*pi*sqrt(a^3/mu_E); % Orbital period [s]

% Match the time span with ephemeris data duration
eph_duration = seconds(utc_time(end) - utc_time(1));
tspan = linspace(0, eph_duration, min(10000, length(utc_time))); % Match resolution

% Set options for the ODE solver
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);

% Perform the integration
[T, Y] = ode113(@(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan, y0, options);

%% Calculate ground track for propagated orbit
fprintf('Calculating propagated ground track...\n');
r = Y(:,1:3);
gtrack_data = zeros(length(T), 4);

for i = 1:length(T)
    current_r = r(i,:);
    R = norm(current_r);
    
    % Calculate declination (latitude in radians)
    delta = asin(current_r(3) / R);
    gtrack_data(i,1) = delta;
    
    % Calculate right ascension
    alpha = atan2(current_r(2), current_r(1));
    gtrack_data(i,2) = alpha;
    
    % Calculate Greenwich sidereal time
    theta_g = theta_g_t0_rad + w_earth_rad * T(i);
    
    % Calculate longitude (wrap to [-180, 180])
    lon = rad2deg(alpha - theta_g);
    gtrack_data(i,3) = mod(lon + 180, 360) - 180;
    
    % Calculate latitude in degrees
    gtrack_data(i,4) = rad2deg(delta);
end

%% Plot comparison
fprintf('Plotting results...\n');

figure(2)
textureImage = imread('EarthTexture.jpg');
imagesc([-180 180], [-90 90], flipud(textureImage));
set(gca, 'YDir', 'normal');
hold on

% Function to plot ground track with proper wrapping
function plot_ground_track(lon, lat, color, style)
    % Detect and handle discontinuities
    wrap_indices = find(abs(diff(lon)) > 180);
    starts = [1; wrap_indices + 1];
    ends = [wrap_indices; length(lon)];
    
    % Plot each continuous segment
    for i = 1:length(starts)
        idx = starts(i):ends(i);
        if length(idx) > 1
            plot(lon(idx), lat(idx), [color style], 'LineWidth', 1.5);
        end
    end
end

% Plot propagated orbit (red line)
plot_ground_track(gtrack_data(:,3), gtrack_data(:,4), 'r', '-');

% Plot ephemeris data (black dots)
plot_ground_track(gtrack_data_eph(:,3), gtrack_data_eph(:,4), 'k', '.');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track Comparison')
legend('Propagated (J2)', 'Ephemeris', 'Location', 'best')
xlim([-180 180])
ylim([-90 90])
grid on
hold off

%% Additional diagnostic plots
figure(3)
subplot(2,1,1)
plot(gtrack_data_eph(:,3), gtrack_data_eph(:,4), 'k.', 'MarkerSize', 8)
hold on
plot(gtrack_data(:,3), gtrack_data(:,4), 'r-', 'LineWidth', 1)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title('Ground Track Comparison')
legend('Ephemeris', 'Propagated (J2)')
grid on

subplot(2,1,2)
% Plot altitude comparison
altitude_eph = zeros(length(utc_time), 1);
for i = 1:length(utc_time)
    [r_eph, ~] = kep2car(keplerian_elements_eph(i,:), mu_E);
    altitude_eph(i) = norm(r_eph) - R_E;
end

altitude_prop = vecnorm(r, 2, 2) - R_E;
time_hours = T/3600; % Convert to hours

plot(time_hours, altitude_prop, 'r-', 'LineWidth', 2)
hold on
plot(seconds(utc_time - utc_time(1))/3600, altitude_eph, 'k.', 'MarkerSize', 8)
xlabel('Time [hours]')
ylabel('Altitude [km]')
title('Altitude Comparison')
legend('Propagated (J2)', 'Ephemeris')
grid on

%% Display statistics
fprintf('\n=== Statistics ===\n');
fprintf('Ephemeris duration: %.2f hours\n', eph_duration/3600);
fprintf('Propagation duration: %.2f hours\n', T(end)/3600);
fprintf('Orbital period: %.2f minutes\n', T_period/60);
fprintf('Number of ephemeris points: %d\n', length(utc_time));
fprintf('Number of propagation points: %d\n', length(T));
