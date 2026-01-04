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
fprintf('Nominal Semi-Major Axis: %4f km\n', a_nominal);
fprintf('Repeating Semi-Major Axis: %4f km\n', a_repeat);
fprintf('Repeating Perturbed Semi-Major Axis: %4f km\n', a_repeat_perturbed);
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
    idx = starts_1orbit(i):ends_1orbit(i);
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

%% Perturbed Orbit for a span of 2 years
% This section integrates the orbit under J2 and drag perturbations for a span
% of 2 years, extracts position/velocity data, and prepares for Keplerian
% elements history computation.
tspan = linspace(0, 2*365*24*60*60, 2*365*24*60); % 2 years in seconds
[~, YPerturbedLong] = ode113(@(t,y) ode_2bp_j2_drag(t,y,mu_E,J2_E, R_E, omega_E, c_d, AreaOverMass), tspan, y_perturbed_initial, options);
r_perturbed_2years = YPerturbedLong(:,1:3);
v_perturbed_2years = YPerturbedLong(:,4:6);
fprintf('Calculating Perturbed Orbit with J2 and Drag for 2 years.\n');
fprintf('----------------------------------------\n');



%% Elements History & Filtering
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
keplerian_history(:,3) = unwrap(keplerian_history(:,3));
keplerian_history(:,4) = unwrap(keplerian_history(:,4));
keplerian_history(:,5) = unwrap(keplerian_history(:,5));
keplerian_history(:,6) = unwrap(keplerian_history(:,6));

% Calculate orbital period and points per orbit
T_orbital = 2*pi*sqrt(kep_parameters(1)^3/mu_E);
points_per_orbit = floor(T_orbital);

% Plot geometric Keplerian elements
figure('Name','Geometric Keplerian Elements History')
subplot(3,1,1)
hold on
plot(tspan./86400, keplerian_history(:,1), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Semi-major Axis [km]')
xlabel('Time [days]')
ylabel('a [km]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan./86400, rad2deg(keplerian_history(:,3)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Inclination [deg]')
xlabel('Time [days]')
ylabel('i [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan./86400, keplerian_history(:,2), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Eccentricity [-]')
xlabel('Time [days]')
ylabel('e [-]')
legend('Location','best')
grid on

% Plot orientation Keplerian elements
figure('Name','Orientation Keplerian Elements History')
subplot(3,1,1)
hold on
plot(tspan./86400, rad2deg(keplerian_history(:,4)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Right Ascension Ascending Node [deg]')
xlabel('Time [days]')
ylabel('Ω [deg]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan./86400, rad2deg(keplerian_history(:,5)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Argument of Periapsis [deg]')
xlabel('Time [days]')
ylabel('ω [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan./86400, rad2deg(keplerian_history(:,6)), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('True Anomaly [deg]')
xlabel('Time [days]')
ylabel('θ [deg]')
legend('Location','best')
grid on

fprintf('Plotting Keplerian Elements History from Perturbed Orbit.\n');
fprintf('----------------------------------------\n');

%% 6.a. Low-Pass Filtering
% This section applies moving average filtering to smooth the Keplerian elements
% history, using window sizes proportional to the orbital period, and plots
% the filtered results alongside the original data.

% Define window sizes for filtering based on orbital period
window_a = floor(points_per_orbit/20);
window_e = floor(points_per_orbit/10);
window_i = floor(points_per_orbit);
window_Omega = floor(points_per_orbit);
window_omega = floor(points_per_orbit);
window_TA = floor(points_per_orbit/50);

% Ensure minimum window sizes
window_a = max(window_a, 1000);
window_e = max(window_e, 2000);
window_i = max(window_i, 5000);
window_Omega = max(window_Omega, 5000);
window_omega = max(window_omega, 5000);
window_TA = max(window_TA, 200);

% Print window sizes
fprintf('Using window sizes:\n');
fprintf('  a: %d points (%.1f min)\n', window_a, window_a/60);
fprintf('  e: %d points (%.1f min)\n', window_e, window_e/60);
fprintf('  i: %d points (%.1f min)\n', window_i, window_i/60);
fprintf('  Ω: %d points (%.1f min)\n', window_Omega, window_Omega/60);
fprintf('  ω: %d points (%.1f min)\n', window_omega, window_omega/60);
fprintf('  θ: %d points (%.1f min)\n', window_TA, window_TA/60);

% Apply moving average filters
a_movmean = movmean(keplerian_history(:,1), window_a, "Endpoints", "fill");
e_movmean = movmean(keplerian_history(:,2), window_e, "Endpoints", "fill");
i_movmean = movmean(keplerian_history(:,3), window_i, "Endpoints", "fill");
Omega_movmean = movmean(keplerian_history(:,4), window_Omega, "Endpoints", "fill");
omega_movmean = movmean(keplerian_history(:,5), window_omega, "Endpoints", "fill");
TA_movmean = movmean(keplerian_history(:,6), window_TA, "Endpoints", "fill");

% Downsample data for plotting if necessary
downsample_factor = floor(length(tspan)/10000);
if downsample_factor > 1
    idx = 1:downsample_factor:length(tspan);
    tspan_plot = tspan(idx);
    a_plot = keplerian_history(idx,1);
    a_movmean_plot = a_movmean(idx);
    e_plot = keplerian_history(idx,2);
    e_movmean_plot = e_movmean(idx);
    i_plot = keplerian_history(idx,3);
    i_movmean_plot = i_movmean(idx);
    Omega_plot = keplerian_history(idx,4);
    Omega_movmean_plot = Omega_movmean(idx);
    omega_plot = keplerian_history(idx,5);
    omega_movmean_plot = omega_movmean(idx);
    TA_plot = keplerian_history(idx,6);
    TA_movmean_plot = TA_movmean(idx);
else
    tspan_plot = tspan;
    a_plot = keplerian_history(:,1);
    a_movmean_plot = a_movmean;
    e_plot = keplerian_history(:,2);
    e_movmean_plot = e_movmean;
    i_plot = keplerian_history(:,3);
    i_movmean_plot = i_movmean;
    Omega_plot = keplerian_history(:,4);
    Omega_movmean_plot = Omega_movmean;
    omega_plot = keplerian_history(:,5);
    omega_movmean_plot = omega_movmean;
    TA_plot = keplerian_history(:,6);
    TA_movmean_plot = TA_movmean;
end

% Plot filtered geometric elements
figure('Name','Geometric Keplerian Elements History')
subplot(3,1,1)
hold on
plot(tspan_plot./86400, a_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, a_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Semi-major Axis [km]')
xlabel('Time [days]')
ylabel('a [km]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot./86400, rad2deg(i_plot), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, rad2deg(i_movmean_plot), 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Inclination [deg]')
xlabel('Time [days]')
ylabel('i [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot./86400, e_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, e_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Eccentricity [-]')
xlabel('Time [days]')
ylabel('e [-]')
legend('Location','best')
grid on

% Plot filtered orientation elements
figure('Name','Orientation Keplerian Elements History - Corrected')
subplot(3,1,1)
hold on
plot(tspan_plot./86400, rad2deg(Omega_plot), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, rad2deg(Omega_movmean_plot), 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Right Ascension Ascending Node [deg]')
xlabel('Time [days]')
ylabel('Ω [deg]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot./86400, rad2deg(omega_plot), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, rad2deg(omega_movmean_plot), 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Argument of Periapsis [deg]')
xlabel('Time [days]')
ylabel('ω [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot./86400, rad2deg(TA_plot), 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot./86400, rad2deg(TA_movmean_plot), 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('True Anomaly [deg]')
xlabel('Time [days]')
ylabel('θ [deg]')
legend('Location','best')
grid on

fprintf('Plotting Filtered Keplerian Elements History from Perturbed Orbit.\n');
fprintf('----------------------------------------\n');

%% 3.b. Integration using Gauss Planetary Equations for a span of 2 years
% This section integrates the Keplerian elements directly using Gauss planetary
% equations with perturbations, providing an alternative to Cartesian integration.

% Initial Keplerian elements as column vector
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
fprintf('----------------------------------------\n');

%% 3.c. Compare with Cartesian integration
% This section compares the results from Gauss and Cartesian integrations by
% calculating position differences and plotting them.

% Calculate differences if time spans match
if length(TPerturbed) == length(TGauss)
    pos_diff = r_perturbed_2years - r_gauss;
    pos_magnitude_diff = sqrt(sum(pos_diff.^2, 2));
    
    figure('Name', 'Comparison: Cartesian vs Gauss Integration');
    subplot(2,1,1);
    plot(TGauss./86400, pos_magnitude_diff);
    xlabel('Time [days]');
    ylabel('Position Difference [km]');
    title('Position Discrepancy between Integration Methods');
    grid on;
    
    subplot(2,1,2);
    semilogy(TGauss./86400, pos_magnitude_diff);
    xlabel('Time [days]');
    ylabel('Position Difference [km] (log scale)');
    title('Logarithmic Scale');
    grid on;
    
    fprintf('Maximum position difference: %.6f km\n', max(pos_magnitude_diff));
    fprintf('Mean position difference: %.6f km\n', mean(pos_magnitude_diff));
end

fprintf('Comparing Gauss and Cartesian Integration Results.\n');
fprintf('----------------------------------------\n');
%% Analyze secular rates from Gauss equations
% This section computes and plots the secular rates of change for Keplerian elements
% from the Gauss integration, filtering them to isolate long-term trends.

% Compute gradients for rates
da_dt = gradient(KepGauss(:,1), TGauss);
de_dt = gradient(KepGauss(:,2), TGauss);
di_dt = gradient(KepGauss(:,3), TGauss);
dOmega_dt = gradient(KepGauss(:,4), TGauss);
domega_dt = gradient(KepGauss(:,5), TGauss);

% Filter rates using moving average
window_secular = floor(points_per_orbit * 10);
da_dt_secular = movmean(da_dt, window_secular, "Endpoints", "fill");
de_dt_secular = movmean(de_dt, window_secular, "Endpoints", "fill");
di_dt_secular = movmean(di_dt, window_secular, "Endpoints", "fill");
dOmega_dt_secular = movmean(dOmega_dt, window_secular, "Endpoints", "fill");
domega_dt_secular = movmean(domega_dt, window_secular, "Endpoints", "fill");

% Plot secular rates
figure('Name', 'Secular Rates from Gauss Equations');
subplot(3,2,1);
plot(TGauss./86400, da_dt_secular * 86400 * 365.25);
xlabel('Time [days]');
ylabel('da/dt [km/year]');
title('Semi-major Axis Rate');
grid on;

subplot(3,2,2);
plot(TGauss./86400, de_dt_secular * 86400 * 365.25);
xlabel('Time [days]');
ylabel('de/dt [1/year]');
title('Eccentricity Rate');
grid on;

subplot(3,2,3);
plot(TGauss./86400, rad2deg(di_dt_secular) * 86400 * 365.25);
xlabel('Time [days]');
ylabel('di/dt [deg/year]');
title('Inclination Rate');
grid on;

subplot(3,2,4);
plot(TGauss./86400, rad2deg(dOmega_dt_secular) * 86400 * 365.25);
xlabel('Time [days]');
ylabel('dΩ/dt [deg/year]');
title('RAAN Rate');
grid on;

subplot(3,2,5);
plot(TGauss./86400, rad2deg(domega_dt_secular) * 86400 * 365.25);
xlabel('Time [days]');
ylabel('dω/dt [deg/year]');
title('Argument of Perigee Rate');
grid on;

fprintf('Analyzing Secular Rates from Gauss Integration.\n');
fprintf('----------------------------------------\n');

%% Theoretical secular rates from J2 (for comparison)
% This section calculates analytical secular rates due to J2 perturbations
% and compares them with numerical results from Gauss integration.

% Extract initial elements
a_initial = kep_parameters(1);
e_initial = kep_parameters(2);
i_initial = kep_parameters(3);
n = sqrt(mu_E / a_initial^3);

% Analytical J2 secular rates
dOmega_dt_J2 = -1.5 * J2_E * (R_E^2) * n * cos(i_initial) / (a_initial^2 * (1 - e_initial^2)^2);
domega_dt_J2 = 0.75 * J2_E * (R_E^2) * n * (5*cos(i_initial)^2 - 1) / (a_initial^2 * (1 - e_initial^2)^2);

fprintf('\nTheoretical secular rates from J2:\n');
fprintf('  dΩ/dt: %.6f deg/day\n', rad2deg(dOmega_dt_J2) * 86400);
fprintf('  dω/dt: %.6f deg/day\n', rad2deg(domega_dt_J2) * 86400);

% Compare with numerical averages
mean_dOmega_num = mean(dOmega_dt_secular);
mean_domega_num = mean(domega_dt_secular);
fprintf('\nNumerical secular rates (average):\n');
fprintf('  dΩ/dt: %.6f deg/day\n', rad2deg(mean_dOmega_num) * 86400);
fprintf('  dω/dt: %.6f deg/day\n', rad2deg(mean_domega_num) * 86400);

fprintf('----------------------------------------\n');



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
subplot(6,1,1)
plot(utc_time_tank,keplerian_elements_eph_tank(:,1))
title('Semi-major Axis')
% 
subplot(6,1,2)
plot(utc_time_tank,keplerian_elements_eph_tank(:,2))
title('Eccentricity')
% 
subplot(6,1,3)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,3)))
title('Inclination')
% 
subplot(6,1,4)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,4)))
title('Right Ascension Ascending Node [deg]')

subplot(6,1,5)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,5)))
title('Argument of Perigee [deg]')
% 
subplot(6,1,6)
plot(utc_time_tank,rad2deg(keplerian_elements_eph_tank(:,6)))
title('True Anomaly [deg]')
% 
% 
% 
figure('Name','Keplerian Elements History for Titan, Ephemeris vs Gaussian Theory')
subplot(6,1,1)
plot(utc_time_titan,keplerian_elements_eph_titan(:,1))
title('Semi-major Axis')
% 
subplot(6,1,2)
plot(utc_time_titan,keplerian_elements_eph_titan(:,2))
title('Eccentricity')
% 
subplot(6,1,3)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,3)))
title('Inclination')
% 
subplot(6,1,4)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,4)))
title('Right Ascension Ascending Node [deg]')
% 
subplot(6,1,5)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,5)))
title('Argument of Perigee [deg]')
% 
subplot(6,1,6)
plot(utc_time_titan,rad2deg(keplerian_elements_eph_titan(:,6)))
title('True Anomaly [deg]')

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
xlim([-9000 9000]); ylim([-9000 9000]); zlim([-9000 9000]);
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
comet_length = 1000; % Length of the comet tail
step = 10; % Step size for faster animation (increase to skip points)
for i = 1:step:length(YPerturbed_1day)
    start_idx = max(1, i - comet_length + 1);
    % Update trajectory
    set(h_traj, 'XData', YPerturbed_1day(start_idx:i,1), 'YData', YPerturbed_1day(start_idx:i,2), 'ZData', YPerturbed_1day(start_idx:i,3));
    
    % Radial arrow: from origin to position
    pos = YPerturbed_1day(i,1:3);
    radial_vec = pos / norm(pos); % Unit vector radial outward
    set(h_radial, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', radial_vec(1)*1000, 'VData', radial_vec(2)*1000, 'WData', radial_vec(3)*1000); % Scaled for visibility
    
    % Perpendicular arrow: normal to velocity in orbital plane (approximate as cross product of position and velocity)
    vel = YPerturbed_1day(i,4:6);
    vel_vevc = vel / norm(vel); % Unit velocity vector
    set(h_vel, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', vel_vevc(1)*1000, 'VData', vel_vevc(2)*1000, 'WData', vel_vevc(3)*1000); % Scaled for visibility

    perp_vec = cross(pos, vel); % Perpendicular to both position and velocity
    perp_vec = perp_vec / norm(perp_vec); % Unit vector
    set(h_perp, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3), 'UData', perp_vec(1)*1000, 'VData', perp_vec(2)*1000, 'WData', perp_vec(3)*1000); % Scaled for visibility
    
    drawnow;
    pause(0.000001); % Short pause; adjust for speed
end

hold off