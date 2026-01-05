clc
clearvars
close all

%% Constants and Setup
% This section initializes the simulation environment by clearing the workspace,
% defining unit conversion functions, adding necessary paths for custom functions
% and resources, and setting up physical constants, Earth model, and initial parameters
% for the orbital mechanics simulation.

% Define conversion functions for degrees to radians and days to seconds
deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;

% Add paths for ephemeris data, custom functions, and time conversion utilities
addpath("Code\Assignment2\Ephemeris\")
addpath("Code\Assignment2\Functions\")
addpath("Code\Assignment2\Functions\timeConversion\")
addpath('Code\Assignment2\img')

% Define physical constants for Earth
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);  % Earth radius [km]
J2_E = astroConstants(33); % Earth J2 coefficient

% Load Earth texture image for visualization
EarthImage = imread('img\earth.jpg');

% Calculate Earth's angular velocity for sidereal rotation
w_earth_rad = 2*pi/(astroConstants(53)*60*60);
omega_E = [0;0;w_earth_rad]; % Angular velocity vector [rad/s]

% Set simulation parameters: number of orbits (m), scaling factor (k), and start date
m = 12;
k = 1;
starting_date = datetime(2037,5,11,13,34,24);

% Compute ground track duration based on Earth's rotation and orbital parameters
T_groundtrack = (2*pi/w_earth_rad)*m/k; % Duration in seconds
ending_date = starting_date + seconds(T_groundtrack);
tspan = starting_date:seconds(1):ending_date; % Time vector for simulation
tspan(end) = []; % Remove last element to avoid duplication

% Define drag parameters: coefficient and area-to-mass ratio
c_d = 2.1; % Drag coefficient [-]
AreaOverMass = 0.01318; % Area-to-mass ratio [m^2/kg]

% Generate a spherical mesh for Earth representation
[XS,YS,ZS] = sphere(50);
scaleFactor = R_E; % Scale factor matching Earth's radius

% Scale the sphere coordinates to Earth's size
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;

% Initialize Earth rotation parameters for ground track calculations
w_earth = 360 / 86164; % Sidereal rotation rate [deg/s]
t_0 = 0; % Reference time
theta_g_t_0 = 0; % Initial Greenwich sidereal time [deg]
theta_g_t0_rad = deg2rad(0); % Initial GST in radians

% Define initial Keplerian orbital elements
kep_parameters = [0.69427e4, 0.0277, deg2rad(85), deg2rad(72), deg2rad(130), 0]; % [a (km), e, i (rad), Omega (rad), omega (rad), theta (rad)]

% Set ODE solver options for high precision
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

fprintf('Setting up. Please wait...\n');
fprintf('----------------------------------------\n');

%% Nominal Orbit
% This section computes and visualizes the unperturbed two-body problem orbit.
% It starts by converting initial Keplerian elements to Cartesian coordinates,
% defines the time span based on the orbital period, integrates the equations
% of motion using ODE113, extracts position and velocity data, and plots the
% orbit in 3D space with Earth as a textured sphere. Finally, it prints key
% orbital characteristics for verification.

% Convert initial Keplerian elements to position and velocity vectors
[r_nominal_initial, v_nominal_initial] = kep2car(kep_parameters, mu_E);
y_nominal_initial = [r_nominal_initial; v_nominal_initial]; % Combine into state vector [r; v]

% Calculate semi-major axis and orbital period for the nominal orbit
a_nominal = kep_parameters(1); % Semi-major axis [km]
T_nominal = 2*pi*sqrt(a_nominal^3/mu_E); % Orbital period [s]

% Define time span for integration: simulate over the ground track duration
% but limit to integer multiples of the period for clarity
tspan_nominal = linspace(0, T_groundtrack, floor(T_groundtrack)); % Time vector [s]

% Integrate the two-body problem equations using ODE113 solver
[TNominal, YNominal] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan_nominal, y_nominal_initial, options);

% Extract position and velocity from the integrated state
r_nominal = YNominal(:,1:3); % Position vectors [km]
v_nominal = YNominal(:,4:6); % Velocity vectors [km/s]

% Create a 3D plot of the nominal orbit with Earth visualization
figure('Name','Nominal Orbit')
hold on
plot3(YNominal(:,1), YNominal(:,2), YNominal(:,3),'r') % Plot the orbit trajectory
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none'); % Earth surface

% Finalize the plot with labels, equal axes, and grid
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbit Representation - Nominal Orbit (2BP)');
axis equal;
grid on;
view(45,30);
legend('Orbit','Earth');

% Print key orbital parameters for the nominal orbit
fprintf('Plotting Nominal Orbit. Identified as a Polar Low Earth Orbit.\n');
fprintf('Radius of the Pericenter: %4f km\n', kep_parameters(1)*(1-kep_parameters(2)));
fprintf('Radius of the Apocenter: %4f km\n', kep_parameters(1)*(1+kep_parameters(2)));
fprintf('Altitude of the Pericenter: %4f km\n', kep_parameters(1)*(1-kep_parameters(2)) - R_E);
fprintf('Altitude of the Apocenter: %4f km\n', kep_parameters(1)*(1+kep_parameters(2)) - R_E);
fprintf('Orbital Period: %4f minutes\n', T_nominal/60);
fprintf('----------------------------------------\n');

%% Ground Track of Nominal Orbit
% This section calculates the ground track for the nominal orbit by converting
% Cartesian positions to geodetic coordinates (longitude and latitude) over time,
% accounting for Earth's rotation. It stores the results in a data array for plotting.

% Initialize array to store ground track data: [delta, alpha, lon, lat]
gtrack_data = zeros(length(TNominal),4);

% Loop through each time step to compute ground track coordinates
for i = 1:length(TNominal)
    R = norm(r_nominal(i,:)); % Magnitude of position vector
    
    % Calculate declination (latitude in radians)
    gtrack_data(i,1) = asin(r_nominal(i,3) / R); % delta
    
    % Calculate right ascension (longitude in inertial frame)
    gtrack_data(i,2) = atan2(r_nominal(i,2), r_nominal(i,1)); % alpha
    
    % Compute Greenwich sidereal time at current time
    theta_g = theta_g_t_0 + w_earth_rad * (TNominal(i) - t_0);
    
    % Calculate longitude by subtracting GST from right ascension, and wrap to [-180, 180]
    lon = rad2deg(gtrack_data(i,2) - theta_g);
    gtrack_data(i,3) = mod(lon + 180, 360) - 180;
    
    % Calculate latitude in degrees
    gtrack_data(i,4) = rad2deg(gtrack_data(i,1));
end

fprintf('Calculating Nominal Ground Track.\n');
fprintf('----------------------------------------\n');
%% Repeating Ground Track
% This section computes a repeating ground track by adjusting the semi-major axis
% to achieve resonance with Earth's rotation, integrates the unperturbed orbit,
% and calculates the corresponding ground track data.

% Calculate semi-major axis for repeating ground track using custom function
a_repeat = RepeatingGroundTrack(m,k,w_earth_rad,mu_E);
kep_repeat = [a_repeat kep_parameters(2) kep_parameters(3) kep_parameters(4) kep_parameters(5) kep_parameters(6)];

% Convert Keplerian elements to initial Cartesian state
[r_repeat_initial, v_repeat_initial] = kep2car(kep_repeat, mu_E);
y_repeat_initial = [r_repeat_initial; v_repeat_initial]; % State vector

% Define time span for repeating orbit integration
tspan_repeat = linspace(0, T_groundtrack, floor(T_groundtrack)); % Time vector

% Integrate the two-body problem for repeating orbit
[TRepeat, YRepeat] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan_repeat, y_repeat_initial, options);

% Extract position and velocity
r_repeat = YRepeat(:,1:3);
v_repeat = YRepeat(:,4:6);

% Initialize ground track data array
gtrack_data_repeat = zeros(length(TNominal),4);

% Compute ground track coordinates for repeating orbit
for i = 1:length(TRepeat)
    R = norm(r_repeat(i,:));
    
    % Calculate declination and right ascension
    gtrack_data_repeat(i,1) = asin(r_repeat(i,3) / R);
    gtrack_data_repeat(i,2) = atan2(r_repeat(i,2), r_repeat(i,1));
    
    % Compute Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TRepeat(i) - t_0);
    
    % Calculate longitude and latitude
    lon = rad2deg(gtrack_data_repeat(i,2) - theta_g);
    gtrack_data_repeat(i,3) = mod(lon + 180, 360) - 180;
    gtrack_data_repeat(i,4) = rad2deg(gtrack_data_repeat(i,1));
end

fprintf('Calculating Repeating Nominal Ground Track assuming k = 1 and m = 12.\n');

fprintf('----------------------------------------\n');

%% Repeating Perturbed Ground Track
% This section computes the ground track for a repeating orbit under perturbations
% (J2 and drag), integrating the perturbed equations and calculating coordinates.

% Use the same semi-major axis as repeating orbit
[a_repeat_perturbed, T_nodal, Omega_dot] = computeRGTSemiMajorAxis(mu_E, R_E, J2_E, w_earth_rad, k, m, kep_parameters(3), kep_parameters(2), 1e-9);
kep_repeat_perturbed = [a_repeat_perturbed kep_parameters(2) kep_parameters(3) kep_parameters(4) kep_parameters(5) kep_parameters(6)];

% Initial state for perturbed repeating orbit
[r_repeat_initial_perturbed, v_repeat_initial_perturbed] = kep2car(kep_repeat_perturbed, mu_E);
y_repeat_initial_perturbed = [r_repeat_initial_perturbed; v_repeat_initial_perturbed];

% Time span for integration
tspan_repeat_perturbed = linspace(0, T_groundtrack, floor(T_groundtrack));

% Integrate with perturbations
[TRepeatPerturbed, YRepeatPerturbed] = ode113(@(t,y) ode_2bp_j2_drag(t,y,mu_E,J2_E, R_E, omega_E, c_d, AreaOverMass), tspan_repeat_perturbed, y_repeat_initial_perturbed, options);

% Extract position and velocity
r_repeat_perturbed = YRepeatPerturbed(:,1:3);
v_repeat_perturbed = YRepeatPerturbed(:,4:6);

% Initialize ground track data
gtrack_data_repeat_perturbed = zeros(length(TRepeatPerturbed),4);

% Compute ground track for perturbed repeating orbit
for i = 1:length(TRepeatPerturbed)
    R = norm(r_repeat_perturbed(i,:));
    
    % Calculate declination and right ascension
    gtrack_data_repeat_perturbed(i,1) = asin(r_repeat_perturbed(i,3) / R);
    gtrack_data_repeat_perturbed(i,2) = atan2(r_repeat_perturbed(i,2), r_repeat_perturbed(i,1));
    
    % Compute Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TRepeatPerturbed(i) - t_0);
    
    % Calculate longitude and latitude
    lon = rad2deg(gtrack_data_repeat_perturbed(i,2) - theta_g);
    gtrack_data_repeat_perturbed(i,3) = mod(lon + 180, 360) - 180;
    gtrack_data_repeat_perturbed(i,4) = rad2deg(gtrack_data_repeat_perturbed(i,1));
end

fprintf('Calculating Repeating Perturbed Ground Track assuming k = 1 and m = 12.\n');
fprintf('Nominal Semi-Major Axis: %.4f km\n', a_nominal);
fprintf('Repeating Semi-Major Axis: %.4f km\n', a_repeat);
fprintf('Repeating Perturbed Semi-Major Axis: %.4f km\n', a_repeat_perturbed);
fprintf('----------------------------------------\n');

%% Perturbed Orbit (J2 + Drag)
% This section integrates the orbit under J2 and atmospheric drag perturbations,
% plots the 3D trajectory, and extracts position/velocity data for further analysis.

% Initial state for perturbed orbit
[r_perturbed_initial, v_perturbed_initial] = kep2car(kep_parameters, mu_E);
y_perturbed_initial = [r_perturbed_initial; v_perturbed_initial];

% Calculate orbital period and time span
a_perturbed_initial = kep_parameters(1);
T_perturbed_initial = 2*pi*sqrt(a_perturbed_initial^3/mu_E);
tspan_perturbed = linspace(0, 5*T_groundtrack, floor(T_groundtrack));

% Integrate with perturbations
[TPerturbed, YPerturbed] = ode113(@(t,y) ode_2bp_j2_drag(t,y,mu_E,J2_E, R_E, omega_E, c_d, AreaOverMass), tspan_perturbed, y_perturbed_initial, options);

% Extract position and velocity
r_perturbed = YPerturbed(:,1:3);
v_perturbed = YPerturbed(:,4:6);

% Plot the perturbed orbit
figure('Name','Perturbed Orbit')
hold on
%plot3(YPerturbed(:,1), YPerturbed(:,2), YPerturbed(:,3))
scatter3(YPerturbed(:,1), YPerturbed(:,2), YPerturbed(:,3), 10, 1:length(YPerturbed), 'filled');
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none');
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
view(45,30);

title('Orbit Representation - Perturbed Orbit (J2 + Drag)');
axis equal;
grid on;
view(45,30);
legend('Orbit','Earth');
colorbar;
caxis([1/(60*24*60) length(YPerturbed)/(60*24*60)]); % Set color axis to represent time in days

% Add colorbar for orbit coloring by index


fprintf('Calculating and Plotting Perturbed Orbit with J2 and Drag.\n');
fprintf('----------------------------------------\n');

%% 2.c. Ground Track of Perturbed
% This section calculates the ground track for the perturbed orbit.

% Initialize ground track data
gtrack_data_perturbed = zeros(length(TNominal),4);

% Compute coordinates
for i = 1:length(TPerturbed)
    R = norm(r_perturbed(i,:));
    
    % Calculate declination and right ascension
    gtrack_data_perturbed(i,1) = asin(r_perturbed(i,3) / R);
    gtrack_data_perturbed(i,2) = atan2(r_perturbed(i,2), r_perturbed(i,1));
    
    % Compute Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TPerturbed(i) - t_0);
    
    % Calculate longitude and latitude
    lon = rad2deg(gtrack_data_perturbed(i,2) - theta_g);
    gtrack_data_perturbed(i,3) = mod(lon + 180, 360) - 180;
    gtrack_data_perturbed(i,4) = rad2deg(gtrack_data_perturbed(i,1));
end

% Handle longitude wrapping for continuous plotting
lon_nominal = gtrack_data(:,3);
lat_nominal = gtrack_data(:,4);
lon_repeat = gtrack_data_repeat(:,3);
lat_repeat = gtrack_data_repeat(:,4);
lon_perturbed = gtrack_data_perturbed(:,3);
lat_perturbed = gtrack_data_perturbed(:,4);
lon_repeat_perturbed = gtrack_data_repeat_perturbed(:,3);
lat_repeat_perturbed = gtrack_data_repeat_perturbed(:,4);

% Detect discontinuities in longitude
wrap_indices_nominal = find(abs(diff(lon_nominal)) > 180);
starts_nominal = [1; wrap_indices_nominal + 1];
ends_nominal = [wrap_indices_nominal; length(lon_nominal)];
wrap_indices_repeat = find(abs(diff(lon_repeat)) > 180);
starts_repeat = [1; wrap_indices_repeat + 1];
ends_repeat = [wrap_indices_repeat; length(lon_repeat)];
wrap_indices_perturbed = find(abs(diff(lon_perturbed)) > 180);
starts_perturbed = [1; wrap_indices_perturbed + 1];
ends_perturbed = [wrap_indices_perturbed; length(lon_perturbed)];
wrap_indices_repeat_perturbed = find(abs(diff(lon_repeat_perturbed)) > 180);
starts_repeat_perturbed = [1; wrap_indices_repeat_perturbed + 1];
ends_repeat_perturbed = [wrap_indices_repeat_perturbed; length(lon_repeat_perturbed)];

fprintf('Plotting Ground Track for Nominal, Repeating, and Perturbed Orbits.\n');
fprintf('----------------------------------------\n');

%% Ground Track Nominal & Perturbed Orbit for 1 period, 1 day, 12 days
% Plot ground track for nominal orbit only - 1 orbit period

% Extract data for 1 orbit period (nominal and perturbed)
idx_1orbit = find(TNominal <= T_nominal);
lon_1orbit = lon_nominal(idx_1orbit);
lat_1orbit = lat_nominal(idx_1orbit);
lon_1orbit_perturbed = lon_perturbed(idx_1orbit);
lat_1orbit_perturbed = lat_perturbed(idx_1orbit);

% Detect discontinuities in longitude for 1 orbit period
wrap_indices_1orbit = find(abs(diff(lon_1orbit)) > 180);
starts_1orbit = [1; wrap_indices_1orbit + 1];
ends_1orbit = [wrap_indices_1orbit; length(lon_1orbit)];
wrap_indices_1orbit_perturbed = find(abs(diff(lon_1orbit_perturbed)) > 180);
starts_1orbit_perturbed = [1; wrap_indices_1orbit_perturbed + 1];
ends_1orbit_perturbed = [wrap_indices_1orbit_perturbed; length(lon_1orbit_perturbed)];

% Extract data for 1 day (86400 seconds)
idx_1day = find(TNominal <= 86400);
lon_1day = lon_nominal(idx_1day);
lat_1day = lat_nominal(idx_1day);
lon_1day_perturbed = lon_perturbed(idx_1day);
lat_1day_perturbed = lat_perturbed(idx_1day);

% Detect discontinuities in longitude for 1 day
wrap_indices_1day = find(abs(diff(lon_1day)) > 180);
starts_1day = [1; wrap_indices_1day + 1];
ends_1day = [wrap_indices_1day; length(lon_1day)];
wrap_indices_1day_perturbed = find(abs(diff(lon_1day_perturbed)) > 180);
starts_1day_perturbed = [1; wrap_indices_1day_perturbed + 1];
ends_1day_perturbed = [wrap_indices_1day_perturbed; length(lon_1day_perturbed)];

% Extract data for 1 orbit period (repeating and repeating perturbed)
idx_1orbit_repeat = find(TNominal <= T_nominal);
lon_1orbit_repeat = lon_repeat(idx_1orbit_repeat);
lat_1orbit_repeat = lat_repeat(idx_1orbit_repeat);
lon_1orbit_repeat_perturbed = lon_repeat_perturbed(idx_1orbit_repeat);
lat_1orbit_repeat_perturbed = lat_repeat_perturbed(idx_1orbit_repeat);

% Detect discontinuities in longitude for 1 orbit period (repeat)
wrap_indices_1orbit_repeat = find(abs(diff(lon_1orbit_repeat)) > 180);
starts_1orbit_repeat = [1; wrap_indices_1orbit_repeat + 1];
ends_1orbit_repeat = [wrap_indices_1orbit_repeat; length(lon_1orbit_repeat)];
wrap_indices_1orbit_repeat_perturbed = find(abs(diff(lon_1orbit_repeat_perturbed)) > 180);
starts_1orbit_repeat_perturbed = [1; wrap_indices_1orbit_repeat_perturbed + 1];
ends_1orbit_repeat_perturbed = [wrap_indices_1orbit_repeat_perturbed; length(lon_1orbit_repeat_perturbed)];

% Extract data for 1 day (repeating and repeating perturbed)
idx_1day_repeat = find(TNominal <= 86400);
lon_1day_repeat = lon_repeat(idx_1day_repeat);
lat_1day_repeat = lat_repeat(idx_1day_repeat);
lon_1day_repeat_perturbed = lon_repeat_perturbed(idx_1day_repeat);
lat_1day_repeat_perturbed = lat_repeat_perturbed(idx_1day_repeat);

% Detect discontinuities in longitude for 1 day (repeat)
wrap_indices_1day_repeat = find(abs(diff(lon_1day_repeat)) > 180);
starts_1day_repeat = [1; wrap_indices_1day_repeat + 1];
ends_1day_repeat = [wrap_indices_1day_repeat; length(lon_1day_repeat)];
wrap_indices_1day_repeat_perturbed = find(abs(diff(lon_1day_repeat_perturbed)) > 180);
starts_1day_repeat_perturbed = [1; wrap_indices_1day_repeat_perturbed + 1];
ends_1day_repeat_perturbed = [wrap_indices_1day_repeat_perturbed; length(lon_1day_repeat_perturbed)];


figure('Name','Ground Track Nominal Orbit - 1 Orbit Period')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_1orbit)
    if i == 1
        plot(lon_1orbit(starts_1orbit(i):ends_1orbit(i)), lat_1orbit(starts_1orbit(i):ends_1orbit(i)), 'r', 'LineWidth', 1, 'DisplayName', 'Nominal Orbit');
    else
        plot(lon_1orbit(starts_1orbit(i):ends_1orbit(i)), lat_1orbit(starts_1orbit(i):ends_1orbit(i)), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_1orbit_perturbed)
    idx = starts_1orbit_perturbed(i):ends_1orbit_perturbed(i);
    if i == 1
        plot(lon_1orbit_perturbed(idx), lat_1orbit_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Perturbed Orbit');
    else
        plot(lon_1orbit_perturbed(idx), lat_1orbit_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_1orbit(1), lat_1orbit(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Nominal Start');
plot(lon_1orbit(end), lat_1orbit(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Nominal End');
plot(lon_1orbit_perturbed(1), lat_1orbit_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Perturbed Start');
plot(lon_1orbit_perturbed(end), lat_1orbit_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Perturbed End');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Nominal Orbit (1 Orbit Period)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

% Plot ground track for nominal orbit only - 1 day
figure('Name','Ground Track Nominal Orbit - 1 Day')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_1day)
    idx = starts_1day(i):ends_1day(i);
    if i == 1
        plot(lon_1day(idx), lat_1day(idx), 'r', 'LineWidth', 1, 'DisplayName', 'Nominal Orbit');
    else
        plot(lon_1day(idx), lat_1day(idx), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_1day_perturbed)
    idx = starts_1day_perturbed(i):ends_1day_perturbed(i);
    if i == 1
        plot(lon_1day_perturbed(idx), lat_1day_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Perturbed Orbit');
    else
        plot(lon_1day_perturbed(idx), lat_1day_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_1day(1), lat_1day(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Nominal Start');
plot(lon_1day(end), lat_1day(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Nominal End');
plot(lon_1day_perturbed(1), lat_1day_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Perturbed Start');
plot(lon_1day_perturbed(end), lat_1day_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Perturbed End');
xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Nominal Orbit (1 Day)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

% Plot ground track for nominal orbit only - 12 days
figure('Name','Ground Track Nominal Orbit - 12 Days')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_nominal)
    idx = starts_nominal(i):ends_nominal(i);
    if i == 1
        plot(lon_nominal(idx), lat_nominal(idx), 'r', 'LineWidth', 1, 'DisplayName', 'Nominal Orbit');
    else
        plot(lon_nominal(idx), lat_nominal(idx), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_perturbed)
    idx = starts_perturbed(i):ends_perturbed(i);
    if i == 1
        plot(lon_perturbed(idx), lat_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Perturbed Orbit');
    else
        plot(lon_perturbed(idx), lat_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_nominal(1), lat_nominal(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Nominal Start');
plot(lon_nominal(end), lat_nominal(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Nominal End');
plot(lon_perturbed(1), lat_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Perturbed Start');
plot(lon_perturbed(end), lat_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Perturbed End');
xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Nominal Orbit (12 Days)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

figure('Name','Ground Track Repeating Orbit - 1 Orbit Period')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_1orbit_repeat)
    idx = starts_1orbit_repeat(i):ends_1orbit_repeat(i);
    if i == 1
        plot(lon_1orbit_repeat(idx), lat_1orbit_repeat(idx), 'r', 'LineWidth', 1, 'DisplayName', 'Repeating Orbit');
    else
        plot(lon_1orbit_repeat(idx), lat_1orbit_repeat(idx), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_1orbit_repeat_perturbed)
    idx = starts_1orbit_repeat_perturbed(i):ends_1orbit_repeat_perturbed(i);
    if i == 1
        plot(lon_1orbit_repeat_perturbed(idx), lat_1orbit_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Repeating Perturbed Orbit');
    else
        plot(lon_1orbit_repeat_perturbed(idx), lat_1orbit_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_1orbit_repeat(1), lat_1orbit_repeat(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Repeating Start');
plot(lon_1orbit_repeat(end), lat_1orbit_repeat(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Repeating End');
plot(lon_1orbit_repeat_perturbed(1), lat_1orbit_repeat_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed Start');
plot(lon_1orbit_repeat_perturbed(end), lat_1orbit_repeat_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed End');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Repeating Orbit (1 Orbit Period)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

figure('Name','Ground Track Repeating Orbit - 1 Day')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_1day_repeat)
    idx = starts_1day_repeat(i):ends_1day_repeat(i);
    if i == 1
        plot(lon_1day_repeat(idx), lat_1day_repeat(idx), 'r', 'LineWidth', 1, 'DisplayName', 'Repeating Orbit');
    else
        plot(lon_1day_repeat(idx), lat_1day_repeat(idx), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_1day_repeat_perturbed)
    idx = starts_1day_repeat_perturbed(i):ends_1day_repeat_perturbed(i);
    if i == 1
        plot(lon_1day_repeat_perturbed(idx), lat_1day_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Repeating Perturbed Orbit');
    else
        plot(lon_1day_repeat_perturbed(idx), lat_1day_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_1day_repeat(1), lat_1day_repeat(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Repeating Start');
plot(lon_1day_repeat(end), lat_1day_repeat(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Repeating End');
plot(lon_1day_repeat_perturbed(1), lat_1day_repeat_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed Start');
plot(lon_1day_repeat_perturbed(end), lat_1day_repeat_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed End');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Repeating Orbit (1 Day)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

figure('Name','Ground Track Repeating Orbit - 12 Days')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on
for i = 1:length(starts_repeat)
    idx = starts_repeat(i):ends_repeat(i);
    if i == 1
        plot(lon_repeat(idx), lat_repeat(idx), 'r', 'LineWidth', 1, 'DisplayName', 'Repeating Orbit');
    else
        plot(lon_repeat(idx), lat_repeat(idx), 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

for i = 1:length(starts_repeat_perturbed)
    idx = starts_repeat_perturbed(i):ends_repeat_perturbed(i);
    if i == 1
        plot(lon_repeat_perturbed(idx), lat_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'DisplayName', 'Repeating Perturbed Orbit');
    else
        plot(lon_repeat_perturbed(idx), lat_repeat_perturbed(idx), 'g', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

plot(lon_repeat(1), lat_repeat(1), 'ro', 'MarkerSize', 12, 'DisplayName', 'Repeating Start');
plot(lon_repeat(end), lat_repeat(end), 'rs', 'MarkerSize', 12, 'DisplayName', 'Repeating End');
plot(lon_repeat_perturbed(1), lat_repeat_perturbed(1), 'go', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed Start');
plot(lon_repeat_perturbed(end), lat_repeat_perturbed(end), 'gs', 'MarkerSize', 12, 'DisplayName', 'Repeating Perturbed End');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track - Repeating Orbit (12 Days)')
xlim([-180 180])
ylim([-90 90])
legend('Location','best')
grid on
hold off

fprintf('Plotting Repeating & Repeating Perturbed Ground Track with different periods.\n');
fprintf('----------------------------------------\n');

fprintf('Plotting Nominal & Perturbed Ground Track with different periods.\n');
fprintf('----------------------------------------\n');

%% Keplerian Elements

% Perturbed Orbit for a span of 2 years
% This section integrates the orbit under J2 and drag perturbations for a span
% of 2 years, extracts position/velocity data, and prepares for Keplerian
% elements history computation.
tspan_2years = linspace(0, 2*365*24*60*60, 2*365*24*60); % 2 years in seconds
tStart_2years = tic;
[TPerturbedLong, YPerturbedLong] = ode113(@(t,y) ode_2bp_j2_drag(t,y,mu_E,J2_E, R_E, omega_E, c_d, AreaOverMass), tspan_2years, y_perturbed_initial, options);
r_perturbed_2years = YPerturbedLong(:,1:3);
v_perturbed_2years = YPerturbedLong(:,4:6);
fprintf('Calculating Perturbed Orbit with J2 and Drag for 2 years using Cartesian.\n');
fprintf('----------------------------------------\n');

% This section computes the history of Keplerian elements from the perturbed orbit
% by converting Cartesian states back to Keplerian at each time step, unwraps the
% true anomaly for continuity, and prepares for plotting and filtering.

% Initialize array for Keplerian elements history
keplerian_history = zeros(length(r_perturbed_2years(:,1)),6);

% Convert each Cartesian state to Keplerian elements
for i=1:length(keplerian_history(:,1))
    kep_temp = car2kep(r_perturbed_2years(i,:),v_perturbed_2years(i,:),mu_E);
    keplerian_history(i,1) = kep_temp(1); % a
    keplerian_history(i,2) = kep_temp(2); % e
    keplerian_history(i,3) = kep_temp(3); % i
    keplerian_history(i,4) = kep_temp(4); % Ω
    keplerian_history(i,5) = kep_temp(5); % ω
    keplerian_history(i,6) = kep_temp(6); % θ
end

% Unwrap true anomaly for smooth plotting
keplerian_history(:,4) = unwrap(keplerian_history(:,4));
keplerian_history(:,5) = unwrap(keplerian_history(:,5));
keplerian_history(:,6) = unwrap(keplerian_history(:,6));
tEnd_2years = toc(tStart_2years);
fprintf('Time taken for 2 years integration using Cartesian: %.2f seconds.\n', tEnd_2years);
fprintf('----------------------------------------\n');
% Calculate orbital period and points per orbit


% Integration using Gauss Planetary Equations for a span of 2 years
% This section integrates the Keplerian elements directly using Gauss planetary
% equations with perturbations, providing an alternative to Cartesian integration.

% Initial Keplerian elements as column vector
tStart_Gauss = tic;
kep_initial = kep_parameters';

% Perturbation function handle
a_per_func = @(t, kep) j2_drag_perturbations_RSW(t, kep, mu_E, J2_E, R_E, omega_E, c_d, AreaOverMass);

% Time span for Gauss integration
tspan_gauss = linspace(0, 2*365*24*60*60, 2*365*24*60); % 2 years in seconds

% Integrate Gauss equations
[TGauss, KepGauss] = ode113(@(t, y) gauss_planetary_equations(t, y, mu_E, a_per_func), tspan_gauss, kep_initial, options);

% Convert Gauss results to Cartesian for comparison
r_gauss = zeros(length(TGauss), 3);
v_gauss = zeros(length(TGauss), 3);
for i = 1:length(TGauss)
    [r_temp, v_temp] = kep2car(KepGauss(i,:)', mu_E);
    r_gauss(i,:) = r_temp';
    v_gauss(i,:) = v_temp';
end

fprintf('Integrating Orbit using Gauss Planetary Equations with J2 and Drag.\n');
tEnd_Gauss = toc(tStart_Gauss);
fprintf('Time taken for 2 years integration using Gauss: %.2f seconds.\n', tEnd_Gauss);
fprintf('----------------------------------------\n');

%% Compare with Cartesian integration
% This section compares the results from Gauss and Cartesian integrations by
% calculating position differences, angular differences, and computing comprehensive
% statistical metrics including RMS error.

% Calculate differences if time spans match
if length(TPerturbedLong) == length(TGauss)
    % Position differences and statistics
    pos_diff = r_perturbed_2years - r_gauss;
    pos_magnitude_diff = sqrt(sum(pos_diff.^2, 2));
    
    % Calculate position RMS error
    pos_rms_error = sqrt(mean(pos_magnitude_diff.^2));
    
    % Absolute differences for semi-major axis and eccentricity (non-angular)
    a_diff_absolute = abs(KepGauss(:,1) - keplerian_history(:,1)); % [km]
    a_diff_relative = a_diff_absolute ./ keplerian_history(:,1) * 100; % [%]
    
    e_diff_absolute = abs(KepGauss(:,2) - keplerian_history(:,2)); % [-]
    e_diff_relative = e_diff_absolute ./ (keplerian_history(:,2) + 1e-6) * 100; % [%] (avoid division by zero)
    
    % Angular differences - calculate minimum angular distance (accounts for cyclicity)
    % For inclination, RAAN, argument of periapsis, and true anomaly
    i_diff_rad = abs(KepGauss(:,3) - keplerian_history(:,3));
    i_diff_rad = min(i_diff_rad, 2*pi - i_diff_rad); % Minimum angular distance
    
    Omega_diff_rad = abs(KepGauss(:,4) - keplerian_history(:,4));
    Omega_diff_rad = min(Omega_diff_rad, 2*pi - Omega_diff_rad);
    
    omega_diff_rad = abs(KepGauss(:,5) - keplerian_history(:,5));
    omega_diff_rad = min(omega_diff_rad, 2*pi - omega_diff_rad);
    
    TA_diff_rad = abs(KepGauss(:,6) - keplerian_history(:,6));
    TA_diff_rad = min(TA_diff_rad, 2*pi - TA_diff_rad);
    
    % Calculate RMS errors for each element
    a_rms_error = sqrt(mean(a_diff_absolute.^2));
    e_rms_error = sqrt(mean(e_diff_absolute.^2));
    i_rms_error = sqrt(mean(i_diff_rad.^2));
    Omega_rms_error = sqrt(mean(Omega_diff_rad.^2));
    omega_rms_error = sqrt(mean(omega_diff_rad.^2));
    TA_rms_error = sqrt(mean(TA_diff_rad.^2));
    
    % Create comprehensive comparison figure
    figure('Name', 'Comparison: Cartesian vs Gauss Integration in Keplerian Elements');
    
    % Semi-major axis
    subplot(2,3,1);
    semilogy(TGauss./86400, a_diff_relative);
    xlabel('Time [days]');
    ylabel('Relative Error [%]');
    title('Semi-major Axis Error');
    grid on;
    
    % Eccentricity
    subplot(2,3,2);
    semilogy(TGauss./86400, e_diff_relative);
    xlabel('Time [days]');
    ylabel('Relative Error [%]');
    title('Eccentricity Error');
    grid on;

    % Inclination
    subplot(2,3,3);
    semilogy(TGauss./86400, rad2deg(i_diff_rad));
    xlabel('Time [days]');
    ylabel('Angular Error [deg]');
    title('Inclination Error');
    grid on;

    % RAAN
    subplot(2,3,4);
    semilogy(TGauss./86400, rad2deg(Omega_diff_rad));
    xlabel('Time [days]');
    ylabel('Angular Error [deg]');
    title('RAAN Error');
    grid on;

    % Argument of periapsis
    subplot(2,3,5);
    semilogy(TGauss./86400, rad2deg(omega_diff_rad));
    xlabel('Time [days]');
    ylabel('Angular Error [deg]');
    title('Argument of Periapsis Error');
    grid on;

    % True Anomaly
    subplot(2,3,6);
    semilogy(TGauss./86400, rad2deg(TA_diff_rad));
    xlabel('Time [days]');
    ylabel('Angular Error [deg]');
    title('True Anomaly Error');
    grid on;
    
    % Print comprehensive statistical analysis
    fprintf('\n========== POSITION STATISTICS ==========\n');
    fprintf('Position Difference [km]:\n');
    fprintf('  Maximum:  %.6e km\n', max(pos_magnitude_diff));
    fprintf('  Minimum:  %.6e km\n', min(pos_magnitude_diff));
    fprintf('  Mean:     %.6e km\n', mean(pos_magnitude_diff));
    fprintf('  Std Dev:  %.6e km\n', std(pos_magnitude_diff));
    fprintf('  RMS:      %.6e km\n', pos_rms_error);
    fprintf('  Median:   %.6e km\n', median(pos_magnitude_diff));
    
    fprintf('\n========== KEPLERIAN ELEMENTS STATISTICS ==========\n');
    
    fprintf('\nSemi-major Axis [km]:\n');
    fprintf('  Absolute Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e\n', ...
        max(a_diff_absolute), min(a_diff_absolute), mean(a_diff_absolute), std(a_diff_absolute), a_rms_error);
    fprintf('  Relative Error %% - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f\n', ...
        max(a_diff_relative), min(a_diff_relative), mean(a_diff_relative), std(a_diff_relative));
    
    fprintf('\nEccentricity [-]:\n');
    fprintf('  Absolute Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e\n', ...
        max(e_diff_absolute), min(e_diff_absolute), mean(e_diff_absolute), std(e_diff_absolute), e_rms_error);
    fprintf('  Relative Error %% - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f\n', ...
        max(e_diff_relative), min(e_diff_relative), mean(e_diff_relative), std(e_diff_relative));
    
    fprintf('\nInclination [rad]:\n');
    fprintf('  Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e [rad]\n', ...
        max(i_diff_rad), min(i_diff_rad), mean(i_diff_rad), std(i_diff_rad), i_rms_error);
    fprintf('  Error - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f [deg]\n', ...
        rad2deg(max(i_diff_rad)), rad2deg(min(i_diff_rad)), rad2deg(mean(i_diff_rad)), rad2deg(std(i_diff_rad)));
    
    fprintf('\nRAAN [rad]:\n');
    fprintf('  Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e [rad]\n', ...
        max(Omega_diff_rad), min(Omega_diff_rad), mean(Omega_diff_rad), std(Omega_diff_rad), Omega_rms_error);
    fprintf('  Error - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f [deg]\n', ...
        rad2deg(max(Omega_diff_rad)), rad2deg(min(Omega_diff_rad)), rad2deg(mean(Omega_diff_rad)), rad2deg(std(Omega_diff_rad)));
    
    fprintf('\nArgument of Periapsis [rad]:\n');
    fprintf('  Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e [rad]\n', ...
        max(omega_diff_rad), min(omega_diff_rad), mean(omega_diff_rad), std(omega_diff_rad), omega_rms_error);
    fprintf('  Error - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f [deg]\n', ...
        rad2deg(max(omega_diff_rad)), rad2deg(min(omega_diff_rad)), rad2deg(mean(omega_diff_rad)), rad2deg(std(omega_diff_rad)));
    
    fprintf('\nTrue Anomaly [rad]:\n');
    fprintf('  Error - Max: %.6e,  Min: %.6e,  Mean: %.6e,  Std: %.6e,  RMS: %.6e [rad]\n', ...
        max(TA_diff_rad), min(TA_diff_rad), mean(TA_diff_rad), std(TA_diff_rad), TA_rms_error);
    fprintf('  Error - Max: %.4f,  Min: %.4f,  Mean: %.4f,  Std: %.4f [deg]\n', ...
        rad2deg(max(TA_diff_rad)), rad2deg(min(TA_diff_rad)), rad2deg(mean(TA_diff_rad)), rad2deg(std(TA_diff_rad)));
    
    fprintf('\n========== ERROR SUMMARY TABLE ==========\n');
    fprintf('Element              | Absolute/Angular Error | RMS Error\n');
    fprintf(repmat('-', 1, 75));
    fprintf('\n');
    fprintf('Semi-major Axis      | %.6e km               | %.6e km\n', mean(a_diff_absolute), a_rms_error);
    fprintf('Eccentricity         | %.6e                   | %.6e\n', mean(e_diff_absolute), e_rms_error);
    fprintf('Inclination          | %.6e rad (%.4f deg)    | %.6e rad\n', mean(i_diff_rad), rad2deg(mean(i_diff_rad)), i_rms_error);
    fprintf('RAAN                 | %.6e rad (%.4f deg)    | %.6e rad\n', mean(Omega_diff_rad), rad2deg(mean(Omega_diff_rad)), Omega_rms_error);
    fprintf('Arg. Periapsis       | %.6e rad (%.4f deg)    | %.6e rad\n', mean(omega_diff_rad), rad2deg(mean(omega_diff_rad)), omega_rms_error);
    fprintf('True Anomaly         | %.6e rad (%.4f deg)    | %.6e rad\n', mean(TA_diff_rad), rad2deg(mean(TA_diff_rad)), TA_rms_error);
    fprintf('Position Vector      | %.6e km                | %.6e km\n', mean(pos_magnitude_diff), pos_rms_error);
end

fprintf('Comparing Gauss and Cartesian Integration Results.\n');
fprintf('----------------------------------------\n');

% Plot geometric Keplerian elements
figure('Name','Geometric Keplerian Elements History')
subplot(2,3,1)
hold on
plot(TPerturbedLong./86400, keplerian_history(:,1), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Semi-major Axis [km]')
xlabel('Time [days]')
ylabel('a [km]')
legend('Location','best')
grid on

subplot(2,3,2)
hold on
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,3)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Inclination [deg]')
xlabel('Time [days]')
ylabel('i [deg]')
legend('Location','best')
grid on

subplot(2,3,3)
hold on
plot(TPerturbedLong./86400, keplerian_history(:,2), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Eccentricity [-]')
xlabel('Time [days]')
ylabel('e [-]')
legend('Location','best')
grid on

subplot(2,3,4)
hold on
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,4)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Right Ascension Ascending Node [deg]')
xlabel('Time [days]')
ylabel('Ω [deg]')
legend('Location','best')
grid on

subplot(2,3,5)
hold on
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,5)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Argument of Periapsis [deg]')
xlabel('Time [days]')
ylabel('ω [deg]')
legend('Location','best')
grid on

subplot(2,3,6)
hold on
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,6)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('True Anomaly [deg]')
xlabel('Time [days]')
ylabel('θ [deg]')
legend('Location','best')
grid on

fprintf('Plotting Keplerian Elements History from Perturbed Orbit (Cartesian Integration).\n');
fprintf('----------------------------------------\n');

%% Low-Pass Filtering
% This section applies moving average filtering to smooth the Keplerian elements
% history, using window sizes determined via Fast Fourier Transform (FFT) analysis.
% Plots FFT spectra, Bode responses, and filtered results.

% Verify data consistency
assert(length(TPerturbedLong) == length(keplerian_history), 'Data length mismatch between TPerturbedLong and keplerian_history');

% Calculate orbital period and points per orbit
T_orbital = 2*pi*sqrt(kep_parameters(1)^3/mu_E);
dt = 60; % Sampling interval in seconds
points_per_orbit = floor(T_orbital / dt);
orbital_freq = 1 / T_orbital; % Orbital frequency in Hz

fprintf('\n=== Low-Pass Filtering Configuration ===\n');
fprintf('Orbital Period: %.2f seconds (%.4f hours)\n', T_orbital, T_orbital/3600);
fprintf('Orbital Frequency: %.6e Hz\n', orbital_freq);
fprintf('Sampling Interval: %d seconds\n', dt);
fprintf('Points per Orbit: %d\n', points_per_orbit);
fprintf('Total data points: %d\n', length(TPerturbedLong));
fprintf('Total duration: %.2f days\n\n', TPerturbedLong(end)/86400);

% Adjust target ratios based on typical perturbation frequencies
% J2 effects: typically 1-10 cycles per day for LEO
% Drag effects: even slower secular decay
% Use 0.01 to 0.1 of orbital frequency for better results
target_ratios = [0.05, 0.05, 0.02, 0.01, 0.01, 0.1]; % [a, e, i, Ω, ω, θ]
element_names = {'Semi-major Axis (a)', 'Eccentricity (e)', 'Inclination (i)', 'RAAN (Ω)', 'A. of Periapsis (ω)', 'True Anomaly (θ)'};

% ========================================
% STEP 1: FFT Analysis and Window Sizing
% ========================================
window_sizes = zeros(1,6);
figure('Name', 'FFT Spectra of Keplerian Elements');

for k = 1:6
    signal = keplerian_history(:,k);
    N = length(signal);
    Fs = 1 / dt; % Sampling frequency in Hz
    
    % Frequency vector up to Nyquist
    f = Fs * (0:(N/2))/N;
    
    % Compute FFT
    Y = fft(signal);
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Calculate cutoff frequency
    cutoff_freq = target_ratios(k) * orbital_freq;
    
    % Calculate window size from cutoff frequency
    window_sizes(k) = max(ceil(1 / (cutoff_freq * dt)), 1);
    
    % Plot FFT spectrum with log scale
    subplot(3,2,k);
    semilogy(f, P1, 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Mark key frequencies
    yrange = ylim;
    plot([orbital_freq, orbital_freq], yrange, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('f_{orb} = %.2e Hz', orbital_freq));
    plot([cutoff_freq, cutoff_freq], yrange, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('f_{cut} = %.2e Hz', cutoff_freq));
    
    hold off;
    set(gca, 'XScale', 'log');
    title(sprintf('%s (Initial window: %d points)', element_names{k}, window_sizes(k)));
    xlabel('Frequency [Hz] (log scale)');
    ylabel('Magnitude (log scale)');
    legend('Spectrum', 'Orbital Freq', 'Cutoff Freq', 'Location', 'best');
    grid on;
end

% ========================================
% STEP 2: Apply Minimum Window Constraints
% ========================================
% Adjust window sizes with physically meaningful minimums
window_a = max(window_sizes(1), 1000);
window_e = max(window_sizes(2), 2000);
window_i = max(window_sizes(3), 5000);
window_Omega = max(window_sizes(4), 5000);
window_omega = max(window_sizes(5), 5000);
window_TA = max(window_sizes(6), 500);

% Validate windows don't exceed data length
max_window = max([window_a, window_e, window_i, window_Omega, window_omega, window_TA]);
assert(max_window < length(keplerian_history), sprintf('Window size %d exceeds data length %d', max_window, length(keplerian_history)));

% Print detailed window information
fprintf('Window Sizes (after minimum constraints):\n');
fprintf('  a (Semi-major axis):    %d points (%.2f hours, ratio: %.4f)\n', window_a, window_a*dt/3600, target_ratios(1));
fprintf('  e (Eccentricity):       %d points (%.2f hours, ratio: %.4f)\n', window_e, window_e*dt/3600, target_ratios(2));
fprintf('  i (Inclination):        %d points (%.2f hours, ratio: %.4f)\n', window_i, window_i*dt/3600, target_ratios(3));
fprintf('  Ω (RAAN):               %d points (%.2f hours, ratio: %.4f)\n', window_Omega, window_Omega*dt/3600, target_ratios(4));
fprintf('  ω (Arg. Periapsis):     %d points (%.2f hours, ratio: %.4f)\n', window_omega, window_omega*dt/3600, target_ratios(5));
fprintf('  θ (True Anomaly):       %d points (%.2f hours, ratio: %.4f)\n', window_TA, window_TA*dt/3600, target_ratios(6));
fprintf('========================================\n\n');

% ========================================
% STEP 3: Bode Plot Analysis
% ========================================
%% Improved Bode Plot Analysis
figure('Name', 'Improved Bode Plots of Moving Average Filters');
Fs = 1 / dt;  % Sampling frequency in Hz

for k = 1:6
    window_sizes_array = [window_a, window_e, window_i, window_Omega, window_omega, window_TA];
    window_size = window_sizes_array(k);
    
    % Create proper filter response using built-in functions
    b = ones(1, window_size) / window_size;  % Moving average coefficients
    a = 1;  % FIR filter
    
    % Compute frequency response properly
    [H, f_bode] = freqz(b, a, 1024, Fs);
    
    % Magnitude in dB
    mag_dB = 20*log10(abs(H) + eps);  % Add eps to avoid log(0)
    
    % Phase in degrees (unwrapped)
    phase_deg = rad2deg(unwrap(angle(H)));
    
    % Find -3dB point
    [~, idx_3db] = min(abs(mag_dB + 3));
    f_3db = f_bode(idx_3db);
    cutoff_freq_k = target_ratios(k) * orbital_freq;
    
    % Magnitude plot
    subplot(6, 2, 2*k-1);
    semilogx(f_bode, mag_dB, 'b-', 'LineWidth', 2);
    hold on;
    
    % Mark key frequencies
    yrange = ylim;
    plot([orbital_freq, orbital_freq], yrange, 'r-', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('f_{orb}=%.2e', orbital_freq));
    plot([cutoff_freq_k, cutoff_freq_k], yrange, 'g-', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('f_{cut}=%.2e', cutoff_freq_k));
    plot(f_3db, -3, 'ko', 'MarkerSize', 10, ...
        'DisplayName', sprintf('f_{-3dB}=%.2e', f_3db));
    
    hold off;
    title(sprintf('Magnitude: %s\n(Window: %d pts)', element_names{k}, window_size));
    ylabel('Magnitude [dB]');
    grid on;
    legend('Location', 'best', 'FontSize', 7);
    set(gca, 'XScale', 'log');
    ylim([-60, 5]);  % Consistent y-axis
    
    % Phase plot
    subplot(6, 2, 2*k);
    semilogx(f_bode, phase_deg, 'r-', 'LineWidth', 2);
    hold on;
    
    % Mark key frequencies
    plot([orbital_freq, orbital_freq], ylim, 'r-', 'LineWidth', 2.5);
    plot([cutoff_freq_k, cutoff_freq_k], ylim, 'g-', 'LineWidth', 2.5);
    plot(f_3db, phase_deg(idx_3db), 'ko', 'MarkerSize', 10);
    
    hold off;
    title(sprintf('Phase: %s', element_names{k}));
    ylabel('Phase [deg]');
    xlabel('Frequency [Hz]');
    grid on;
    set(gca, 'XScale', 'log');
    
    % Store diagnostics
    bode_diagnostics{k} = struct(...
        'f_orbital', orbital_freq, ...
        'f_cutoff', cutoff_freq_k, ...
        'window_size', window_size);
end

% ========================================
% PRINT BODE DIAGNOSTICS
% ========================================
fprintf('Bode Plot Diagnostics (Method 3):\n');
fprintf('%-20s | f_orbital  | f_cutoff   | Window\n', 'Element');
fprintf('%-20s | (Hz)       | (Hz)       | (pts)\n', '');
fprintf('----------------------------------+------------+------------+-------\n');

for k = 1:6
    diag = bode_diagnostics{k};
    fprintf('%-20s | %.2e | %.2e | %d\n', ...
        element_names{k}, diag.f_orbital, diag.f_cutoff, diag.window_size);
end

fprintf('\n✓ Ratio f_-3dB/f_cutoff should be close to 1.0 for proper filter design\n');
fprintf('✓ Method 3 automatically detected cutoff frequencies from spectral slope changes\n');
fprintf('========================================\n\n');

% ========================================
% STEP 4: Apply Filters
% ========================================
fprintf('Applying moving average filters...\n');

a_movmean = movmean(keplerian_history(:,1), window_a, "Endpoints", "fill");
e_movmean = movmean(keplerian_history(:,2), window_e, "Endpoints", "fill");
i_movmean = movmean(keplerian_history(:,3), window_i, "Endpoints", "fill");
Omega_movmean = movmean(keplerian_history(:,4), window_Omega, "Endpoints", "fill");
omega_movmean = movmean(keplerian_history(:,5), window_omega, "Endpoints", "fill");
TA_movmean = movmean(keplerian_history(:,6), window_TA, "Endpoints", "fill");

fprintf('✓ Filters applied successfully\n');

% ========================================
% STEP 5: Plot Filtered Results
% ========================================


% Plot filtered results
figure('Name','Filtered Keplerian Elements History 1');

plot(TPerturbedLong./86400, keplerian_history(:,1), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, a_movmean, 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('Semi-major Axis');
xlabel('Time [days]');
ylabel('[km]');
legend('Location','best');
grid on;

figure('Name','Filtered Keplerian Elements History 2');
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,3)), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, rad2deg(i_movmean), 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('Inclination');
xlabel('Time [days]');
ylabel('[deg]');
legend('Location','best');
grid on;

figure('Name','Filtered Keplerian Elements History 3');
plot(TPerturbedLong./86400, keplerian_history(:,2), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, e_movmean, 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('Eccentricity');
xlabel('Time [days]');
ylabel('[-]');
legend('Location','best');
grid on;

figure('Name','Filtered Keplerian Elements History 4');
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,4)), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, rad2deg(Omega_movmean), 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('Right Ascension Ascending Node');
xlabel('Time [days]');
ylabel('[deg]');
legend('Location','best');
grid on;

figure('Name','Filtered Keplerian Elements History 5');
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,5)), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, rad2deg(omega_movmean), 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('Argument of Periapsis');
xlabel('Time [days]');
ylabel('[deg]');
legend('Location','best');
grid on;

figure('Name','Filtered Keplerian Elements History 6');
plot(TPerturbedLong./86400, rad2deg(keplerian_history(:,6)), 'b', 'DisplayName','Raw', 'LineWidth', 0.5);
hold on;
plot(TPerturbedLong./86400, rad2deg(TA_movmean), 'r', 'DisplayName','Filtered', 'LineWidth', 2);
hold off;
title('True Anomaly');
xlabel('Time [days]');
ylabel('[deg]');
legend('Location','best');
grid on;

fprintf('\n✓ Low-Pass Filtering completed successfully!\n');
fprintf('========================================\n\n');

%% Evolution of two real satellites
% This section analyzes the orbital evolution of two real satellites (BREEZE-M TANK and TITAN 3C TRANSTAGE DEB)
% over a one-month period by reading ephemeris data, propagating with perturbations using Gauss equations,
% computing ground tracks, and plotting comparisons of Keplerian elements and ground tracks.

% Read ephemeris files for the two satellites
[utc_time_tank, mission_info_tank] = read_ephemeris_file('Ephemeris_BreezeM_Tank.txt');
[utc_time_titan, mission_info_titan] = read_ephemeris_file('Ephemeris_Titan3C_Transtage.txt');

% Extract Keplerian elements from ephemeris data
a_eph_tank = mission_info_tank(:,10);
e_eph_tank = mission_info_tank(:,1);
i_eph_tank = deg2rad(mission_info_tank(:,3));
OM_eph_tank = deg2rad(mission_info_tank(:,4));
w_eph_tank = deg2rad(mission_info_tank(:,5));
TA_eph_tank = deg2rad(mission_info_tank(:,9));
keplerian_elements_eph_tank = [a_eph_tank, e_eph_tank, i_eph_tank, OM_eph_tank, w_eph_tank, TA_eph_tank];
keplerian_elements_eph_tank(:,6) = unwrap(keplerian_elements_eph_tank(:,6));

a_eph_titan = mission_info_titan(:,10);
e_eph_titan = mission_info_titan(:,1);
i_eph_titan = deg2rad(mission_info_titan(:,3));
OM_eph_titan = deg2rad(mission_info_titan(:,4));
w_eph_titan = deg2rad(mission_info_titan(:,5));
TA_eph_titan = deg2rad(mission_info_titan(:,9));
keplerian_elements_eph_titan = [a_eph_titan, e_eph_titan, i_eph_titan, OM_eph_titan, w_eph_titan, TA_eph_titan];
keplerian_elements_eph_titan(:,6) = unwrap(keplerian_elements_eph_titan(:,6));

% Initial elements for propagation
kep_initial_tank = keplerian_elements_eph_tank(1,:);
kep_initial_titan = keplerian_elements_eph_titan(1,:);

% Time spans for propagation
tspan_eph_tank = seconds(utc_time_tank - utc_time_tank(1));
tspan_eph_titan = seconds(utc_time_titan - utc_time_titan(1));

% Propagate using Gauss equations
[TGauss_tank, KepGauss_tank] = ode113(@(t, y) gauss_planetary_equations(t, y, mu_E, a_per_func), tspan_eph_tank, kep_initial_tank', options);
[TGauss_titan, KepGauss_titan] = ode113(@(t, y) gauss_planetary_equations(t, y, mu_E, a_per_func), tspan_eph_titan, kep_initial_titan', options);

% Compute ground tracks from Gauss results
gtrack_data_gauss_tank = zeros(length(TGauss_tank), 4);
gtrack_data_gauss_titan = zeros(length(TGauss_titan), 4);
for i = 1:length(TGauss_tank)
    [r_gauss_tank, ~] = kep2car(KepGauss_tank(i,:)', mu_E);
    R_gauss_tank = norm(r_gauss_tank);
    delta_gauss_tank = asin(r_gauss_tank(3) / R_gauss_tank);
    alpha_gauss_tank = atan2(r_gauss_tank(2), r_gauss_tank(1));
    theta_g_gauss_tank = theta_g_t0_rad + w_earth_rad * TGauss_tank(i);
    lon_gauss_tank = rad2deg(alpha_gauss_tank - theta_g_gauss_tank);
    gtrack_data_gauss_tank(i,3) = mod(lon_gauss_tank + 180, 360) - 180;
    gtrack_data_gauss_tank(i,4) = rad2deg(delta_gauss_tank);
end
for i = 1:length(TGauss_titan)
    [r_gauss_titan, ~] = kep2car(KepGauss_titan(i,:)', mu_E);
    R_gauss_titan = norm(r_gauss_titan);
    delta_gauss_titan = asin(r_gauss_titan(3) / R_gauss_titan);
    alpha_gauss_titan = atan2(r_gauss_titan(2), r_gauss_titan(1));
    theta_g_gauss_titan = theta_g_t0_rad + w_earth_rad * TGauss_titan(i);
    lon_gauss_titan = rad2deg(alpha_gauss_titan - theta_g_gauss_titan);
    gtrack_data_gauss_titan(i,3) = mod(lon_gauss_titan + 180, 360) - 180;
    gtrack_data_gauss_titan(i,4) = rad2deg(delta_gauss_titan);
end

fprintf('Processing ephemeris ground track for Tank...\n');
gtrack_data_eph_tank = zeros(length(utc_time_tank), 4);
for i = 1:length(utc_time_tank)
    [r_eph_tank, v_eph_tank] = kep2car(keplerian_elements_eph_tank(i,:), mu_E);
    dt_tank = seconds(utc_time_tank(i) - utc_time_tank(1));
    R_eph_tank = norm(r_eph_tank);
    delta_eph_tank = asin(r_eph_tank(3) / R_eph_tank);
    gtrack_data_eph_tank(i,1) = delta_eph_tank;
    alpha_eph_tank = atan2(r_eph_tank(2), r_eph_tank(1));
    gtrack_data_eph_tank(i,2) = alpha_eph_tank;
    theta_g_eph_tank = theta_g_t0_rad + w_earth_rad * dt_tank;
    lon_eph_tank = rad2deg(alpha_eph_tank - theta_g_eph_tank);
    gtrack_data_eph_tank(i,3) = mod(lon_eph_tank + 180, 360) - 180;
    gtrack_data_eph_tank(i,4) = rad2deg(delta_eph_tank);
end


% Process ephemeris ground track
fprintf('Processing ephemeris ground track for Titan...\n');
gtrack_data_eph_titan = zeros(length(utc_time_titan), 4);
for i = 1:length(utc_time_titan)
    [r_eph_titan, v_eph_titan] = kep2car(keplerian_elements_eph_titan(i,:), mu_E);
    dt_titan = seconds(utc_time_titan(i) - utc_time_titan(1));
    R_eph_titan = norm(r_eph_titan);
    delta_eph_titan = asin(r_eph_titan(3) / R_eph_titan);
    gtrack_data_eph_titan(i,1) = delta_eph_titan;
    alpha_eph_titan = atan2(r_eph_titan(2), r_eph_titan(1));
    gtrack_data_eph_titan(i,2) = alpha_eph_titan;
    theta_g_eph_titan = theta_g_t0_rad + w_earth_rad * dt_titan;
    lon_eph_titan = rad2deg(alpha_eph_titan - theta_g_eph_titan);
    gtrack_data_eph_titan(i,3) = mod(lon_eph_titan + 180, 360) - 180;
    gtrack_data_eph_titan(i,4) = rad2deg(delta_eph_titan);
end

fprintf('Plotting results...\n');

figure('Name','Keplerian Elements History for Tank, Ephemeris vs Gaussian Theory')
subplot(2,3,1)
plot(utc_time_tank,keplerian_elements_eph_tank(:,1))
title('Semi-major Axis [km]')
xlabel('Time')
ylabel('a [km]')
grid on
% 
subplot(2,3,2)
plot(utc_time_tank,keplerian_elements_eph_tank(:,2))
title('Eccentricity [-]')
xlabel('Time')
ylabel('e [-]')
grid on
% 
subplot(2,3,3)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,3)))
title('Inclination [deg]')
xlabel('Time')
ylabel('i [deg]')
grid on
% 
subplot(2,3,4)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,4)))
title('Right Ascension Ascending Node [deg]')
xlabel('Time')
ylabel('Ω [deg]')
%

subplot(2,3,5)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,5)))
title('Argument of Perigee [deg]')
xlabel('Time')
ylabel('ω [deg]')
grid on
% 
subplot(2,3,6)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,6)))
title('True Anomaly [deg]')
xlabel('Time')
ylabel('θ [deg]')
grid on
% 
% 
% 
figure('Name','Keplerian Elements History for Titan, Ephemeris vs Gaussian Theory')
subplot(2,3,1)
plot(utc_time_titan,keplerian_elements_eph_titan(:,1))
title('Semi-major Axis [km]')
xlabel('Time')
ylabel('a [km]')
grid on
% 
subplot(2,3,2)
plot(utc_time_titan,keplerian_elements_eph_titan(:,2))
title('Eccentricity [-]')
xlabel('Time')
ylabel('e [-]')
grid on
% 
subplot(2,3,3)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,3)))
title('Inclination [deg]')
xlabel('Time')
ylabel('i [deg]')
grid on
% 
subplot(2,3,4)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,4)))
title('Right Ascension Ascending Node [deg]')
xlabel('Time')
ylabel('Ω [deg]')
grid on
% 
subplot(2,3,5)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,5)))
title('Argument of Perigee [deg]')
xlabel('Time')
ylabel('ω [deg]')
grid on
% 
subplot(2,3,6)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,6)))
title('True Anomaly [deg]')
xlabel('Time')
ylabel('θ [deg]')
grid on

%% Animation
% This section creates an animated visualization of the perturbed orbit trajectory
% over the Earth sphere for 1 day, using a loop to plot points incrementally for controlled speed.
% Added radial arrow (from Earth center to satellite) and perpendicular arrow (normal to velocity in orbital plane).

% Select data for 1 day (86400 seconds)
idx_1day_anim = find(TPerturbed <= 86400);
YPerturbed_1day = YPerturbed(idx_1day_anim, :);
TPerturbed_1day = TPerturbed(idx_1day_anim);

figure('Name','Animation - 1 Day')
hold on
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
xlim([-9000 9000]);
ylim([-9000 9000]);
zlim([-9000 9000]);

title('Orbit Animation - 1 Day with J2 and Drag Perturbations');
axis equal;
grid on;
view(45,30);
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none');

% Custom animation loop for faster speed control
h_traj = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 2); % Trajectory handle
h_radial = quiver3(NaN, NaN, NaN, NaN, NaN, NaN, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Radial arrow handle
h_vel = quiver3(NaN, NaN, NaN, NaN, NaN, NaN, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Velocity arrow handle
h_perp = quiver3(NaN, NaN, NaN, NaN, NaN, NaN, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Perpendicular arrow handle
comet_length = 300; % Length of the comet tail
step = 10; % Step size for faster animation (increase to skip points)
for i = 1:step:length(YPerturbed_1day)
    start_idx = max(1, i - comet_length + 1);
    % Update trajectory
    set(h_traj, 'XData', YPerturbed_1day(start_idx:i,1), 'YData', YPerturbed_1day(start_idx:i,2), 'ZData', YPerturbed_1day(start_idx:i,3));
    
    % Radial arrow: from origin to position
    pos = YPerturbed_1day(i,1:3);
    vel = YPerturbed_1day(i,4:6);
    radial_vec = pos / norm(pos); % Unit vector radial outward
    set(h_radial, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', radial_vec(1)*1000, 'VData', radial_vec(2)*1000, 'WData', radial_vec(3)*1000); % Scaled for visibility
    
    perp_vec = cross(pos, vel); % Perpendicular to both position and velocity
    perp_vec = perp_vec / norm(perp_vec); % Unit vector
    set(h_perp, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', perp_vec(1)*1000, 'VData', perp_vec(2)*1000, 'WData', perp_vec(3)*1000); % Scaled for visibility

    % Perpendicular arrow: normal to velocity in orbital plane (approximate as cross product of position and velocity)
    vel_vec = cross(perp_vec, pos); % Approximate velocity direction
    vel_vec = vel_vec / norm(vel_vec); % Unit velocity vector
    set(h_vel, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', vel_vec(1)*1000, 'VData', vel_vec(2)*1000, 'WData', vel_vec(3)*1000); % Scaled for visibility

    
    drawnow;
    pause(0.000001); % Short pause; adjust for speed
end

hold off