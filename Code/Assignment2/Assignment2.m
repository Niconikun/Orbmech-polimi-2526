clc
clearvars
close all

%% Constants and Setup
deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;


% Add paths
addpath("Ephemeris\")
addpath("Functions\")
addpath("Functions\timeConversion\")

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);  % Earth radius [km]
J2_E = astroConstants(33); % Earth J2

EarthImage = imread('img\earth.jpg');

w_earth_rad = 2*pi/(astroConstants(53)*60*60);

omega_E = [0;0;w_earth_rad]; %rad/s for a sidereal day

m = 12;
k = 1;
starting_date = datetime(2037,5,11,13,34,24);

T_groundtrack = (2*pi/w_earth_rad)*m/k; % in seconds
ending_date = starting_date + seconds(T_groundtrack);
tspan = starting_date:seconds(1):ending_date;
tspan(end) = [];

c_d = 2.1; %[-] Drag
AreaOverMass = 0.01318; %[m^2/kg]
% Create the sphere
[XS,YS,ZS] = sphere(50);

scaleFactor = R_E; % Radius of Earth in kilometers; adjust as needed

% Scale the sphere
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;


% Earth rotation parameters (sidereal rotation)
w_earth = 360 / 86164; % deg/s (sidereal rotation rate)
t_0 = 0;
theta_g_t_0 = 0; % Initial Greenwich sidereal time [deg]
theta_g_t0_rad = deg2rad(0);    % initial GST in rad

% Initial condition
kep_parameters = [0.69427e4, 0.0277, deg2rad(85), deg2rad(72), deg2rad(130), 0]; % a [km],e [-],i [deg],Omega [deg],omega [deg],theta [deg]


%% 1. Nominal Orbit

[r_nominal_initial, v_nominal_initial] = kep2car(kep_parameters, mu_E);
y_nominal_initial = [r_nominal_initial; v_nominal_initial];  % Combine position and velocity into 6-element vector

% Set time span
a_nominal = kep_parameters(1); % Semi-major axis [km]
T_nominal = 2*pi*sqrt(a_nominal^3/mu_E); % Orbital period [s]
tspan_nominal = linspace( 0, T_groundtrack, floor(T_groundtrack)); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ TNominal, YNominal ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_nominal, y_nominal_initial, options );

% Calculate ground track
r_nominal = YNominal(:,1:3);
v_nominal = YNominal(:,4:6);

% Plot the results
figure('Name','Nominal Orbit')

hold on
plot3( YNominal(:,1), YNominal(:,2), YNominal(:,3))
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none');


% Finalize the plot
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;



%% 2.a. Ground Track of Nominal Orbit

gtrack_data = zeros(length(TNominal),4);

% In your main script, before the ground track loop:

for i = 1:length(TNominal)
    R = norm(r_nominal(i,:));
    %fprintf('Position magnitude: %f km\n', R); % Debug
    
    % Calculate declination (latitude)
    gtrack_data(i,1) = asin(r_nominal(i,3) / R); %delta
    delta = gtrack_data(i,1);
    %fprintf('Declination (delta): %f rad, %f deg\n', delta, rad2deg(delta));
    
    % Calculate right ascension
    gtrack_data(i,2) = atan2(r_nominal(i,2), r_nominal(i,1)); %alpha
    alpha = gtrack_data(i,2);
    %fprintf('Right ascension (alpha): %f rad, %f deg\n', alpha, rad2deg(alpha));
    
    % Calculate Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TNominal(i) - t_0);
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

%% 2.b. Repeating Ground Track

a_repeat = RepeatingGroundTrack(m,k,w_planet_rad,mu_E);
kep_repeat = [a_repeat kep_parameters(2) kep_parameters(3) kep_parameters(4) kep_parameters(5) kep_parameters(6)];

[r_repeat_initial, v_repeat_initial] = kep2car(kep_repeat, mu_E);
y_repeat_initial = [r_repeat_initial; v_repeat_initial];  % Combine position and velocity into 6-element vector

% Set time span

tspan_repeat = linspace( 0, T_groundtrack, floor(T_groundtrack)); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ TRepeat, YRepeat ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_repeat, y_repeat_initial, options );

% Calculate ground track
r_repeat = YRepeat(:,1:3);
v_repeat = YRepeat(:,4:6);

gtrack_data_repeat = zeros(length(TNominal),4);

% In your main script, before the ground track loop:

for i = 1:length(TRepeat)
    R = norm(r_repeat(i,:));
    %fprintf('Position magnitude: %f km\n', R); % Debug
    
    % Calculate declination (latitude)
    gtrack_data_repeat(i,1) = asin(r_repeat(i,3) / R); %delta
    delta = gtrack_data_repeat(i,1);
    %fprintf('Declination (delta): %f rad, %f deg\n', delta, rad2deg(delta));
    
    % Calculate right ascension
    gtrack_data_repeat(i,2) = atan2(r_repeat(i,2), r_repeat(i,1)); %alpha
    alpha = gtrack_data_repeat(i,2);
    %fprintf('Right ascension (alpha): %f rad, %f deg\n', alpha, rad2deg(alpha));
    
    % Calculate Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TRepeat(i) - t_0);
    %fprintf('Greenwich sidereal time (theta_g): %f rad, %f deg\n', theta_g, rad2deg(theta_g));
    
    % Calculate longitude
    lon = rad2deg(alpha - theta_g);
    gtrack_data_repeat(i,3) = mod(lon + 180, 360) - 180;
    lon = gtrack_data_repeat(i,3);
    
    % Calculate latitude
    gtrack_data_repeat(i,4) = rad2deg(delta);
    lat = gtrack_data_repeat(i,4);
    
    %fprintf('Final output - Lon: %f deg, Lat: %f deg\n\n', lon, lat);
end



%% 3.a. Perturbed Orbit (J2 + Drag)

[r_perturbed_initial, v_perturbed_initial] = kep2car(kep_parameters, mu_E);
y_perturbed_initial = [r_perturbed_initial; v_perturbed_initial];  % Combine position and velocity into 6-element vector

% Set time span
a_perturbed_initial = kep_parameters(1); % Semi-major axis [km]
T_perturbed_initial = 2*pi*sqrt(a_perturbed_initial^3/mu_E); % Orbital period [s]
tspan_perturbed = linspace( 0, T_groundtrack, floor(T_groundtrack) ); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ TPerturbed, YPerturbed ] = ode113( @(t,y) ode_2bp_j2_drag(t,y,mu_E,J2_E, R_E, omega_E, c_d, AreaOverMass), tspan_perturbed, y_perturbed_initial, options );

% Calculate ground track
r_perturbed = YPerturbed(:,1:3);
v_perturbed = YPerturbed(:,4:6);

% Plot the results
figure('Name','Perturbed Orbit')

hold on
plot3( YPerturbed(:,1), YPerturbed(:,2), YPerturbed(:,3))
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none');


% Finalize the plot
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;


%% 2.c. Ground Track of Perturbed

gtrack_data_perturbed = zeros(length(TNominal),4);

% In your main script, before the ground track loop:

for i = 1:length(TPerturbed)
    R = norm(r_perturbed(i,:));
    %fprintf('Position magnitude: %f km\n', R); % Debug
    
    % Calculate declination (latitude)
    gtrack_data_perturbed(i,1) = asin(r_perturbed(i,3) / R); %delta
    delta = gtrack_data_perturbed(i,1);
    %fprintf('Declination (delta): %f rad, %f deg\n', delta, rad2deg(delta));
    
    % Calculate right ascension
    gtrack_data_perturbed(i,2) = atan2(r_perturbed(i,2), r_perturbed(i,1)); %alpha
    alpha = gtrack_data_pertubed(i,2);
    %fprintf('Right ascension (alpha): %f rad, %f deg\n', alpha, rad2deg(alpha));
    
    % Calculate Greenwich sidereal time
    theta_g = theta_g_t_0 + w_earth_rad * (TPerturbed(i) - t_0);
    %fprintf('Greenwich sidereal time (theta_g): %f rad, %f deg\n', theta_g, rad2deg(theta_g));
    
    % Calculate longitude
    lon = rad2deg(alpha - theta_g);
    gtrack_data_perturbed(i,3) = mod(lon + 180, 360) - 180;
    lon = gtrack_data_perturbed(i,3);
    
    % Calculate latitude
    gtrack_data_perturbed(i,4) = rad2deg(delta);
    lat = gtrack_data_perturbed(i,4);
    
    %fprintf('Final output - Lon: %f deg, Lat: %f deg\n\n', lon, lat);
end


figure('Name','Ground Track 2BP & Repeating w/ Perturbations')
imagesc([-180 180], [-90 90], flipud(EarthImage));
set(gca, 'YDir', 'normal');
hold on

% Handle longitude wrapping for continuous lines
lon_nominal = gtrack_data(:,3);
lat_nominal = gtrack_data(:,4);

lon_repeat = gtrack_data_repeat(:,3);
lat_repeat = gtrack_data_repeat(:,4);

lon_perturbed = gtrack_data_perturbed(:,3);
lat_perturbed = gtrack_data_perturbed(:,4);


% Detect and handle discontinuities
wrap_indices_nominal = find(abs(diff(lon_nominal)) > 180);
starts_nominal = [1; wrap_indices_nominal + 1];
ends_nominal = [wrap_indices_nominal; length(lon_nominal)];

wrap_indices_repeat = find(abs(diff(lon_repeat)) > 180);
starts_repeat = [1; wrap_indices_repeat + 1];
ends_repeat = [wrap_indices_repeat; length(lon_repeat)];

wrap_indices_perturbed = find(abs(diff(lon_perturbed)) > 180);
starts_perturbed = [1; wrap_indices_perturbed + 1];
ends_perturbed = [wrap_indices_perturbed; length(lon_perturbed)];

% Plot each continuous segment
for i = 1:length(starts_nominal)
    idx = starts_nominal(i):ends(i);
    plot(lon_nominal(idx), lat_nominal(idx), 'r', 'LineWidth', 1);
end

% Plot each continuous segment
for i = 1:length(starts_repeat)
    idx = starts_repeat(i):ends_repeat(i);
    plot(lon_repeat(idx), lat_repeat(idx), 'b', 'LineWidth', 1.5);
end

% Plot each continuous segment
for i = 1:length(starts_perturbed)
    idx = starts_perturbed(i):ends_perturbed(i);
    plot(lon_perturbed(idx), lat_perturbed(idx), 'b', 'LineWidth', 1.5);
end

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track')
xlim([-180 180])
ylim([-90 90])
grid on
hold off


%% 4.a. Elements History & Filtering

keplerian_history = zeros(length(r_perturbed(:,1)),6);

for i=1:length(keplerian_history(:,1))
    kep_temp = car2kep(r_perturbed(i,:),v_perturbed(i,:),mu_E);
    keplerian_history(i,1) = kep_temp(1); % a
    keplerian_history(i,2) = kep_temp(2); % e
    keplerian_history(i,3) = kep_temp(3); % i
    keplerian_history(i,4) = kep_temp(4); % Ω
    keplerian_history(i,5) = kep_temp(5); % ω
    keplerian_history(i,6) = kep_temp(6); % θ
end

keplerian_history(:,6) = unwrap(keplerian_history(:,6));

% Calculate orbital period from initial conditions
T_orbital = 2*pi*sqrt(kep_parameters(1)^3/mu_E); % [s]

% Calculate approximate points per orbit
points_per_orbit = floor(T_orbital); % Since timestep = 1s

% Now create your plots using the downsampled data if applicable
figure('Name','Geometric Keplerian Elements History - Corrected')
subplot(3,1,1)
hold on
plot(tspan_plot, a_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Semi-major Axis [km]')
xlabel('Time')
ylabel('a [km]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot, i_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Inclination [deg]')
xlabel('Time')
ylabel('i [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot, e_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Eccentricity [-]')
xlabel('Time')
ylabel('e [-]')
legend('Location','best')
grid on

figure('Name','Orientation Keplerian Elements History - Corrected')
subplot(3,1,1)
hold on
plot(tspan_plot, Omega_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Right Ascension Ascending Node [deg]')
xlabel('Time')
ylabel('Ω [deg]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot, omega_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('Argument of Periapsis [deg]')
xlabel('Time')
ylabel('ω [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot, TA_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
hold off
title('True Anomaly [deg]')
xlabel('Time')
ylabel('θ [deg]')
legend('Location','best')
grid on

%% 6.a. Low-Pass Filtering

% Use window sizes proportional to orbital period:
window_a = floor(points_per_orbit/20);  % ~5% of orbit for semi-major axis
window_e = floor(points_per_orbit/10);  % ~10% of orbit for eccentricity
window_i = floor(points_per_orbit);     % Full orbit for inclination
window_Omega = floor(points_per_orbit); % Full orbit for RAAN
window_omega = floor(points_per_orbit); % Full orbit for argument of perigee
window_TA = floor(points_per_orbit/50); % Small window for true anomaly

% Ensure windows are reasonable (not too small for large dataset)
window_a = max(window_a, 1000);
window_e = max(window_e, 2000);
window_i = max(window_i, 5000);
window_Omega = max(window_Omega, 5000);
window_omega = max(window_omega, 5000);
window_TA = max(window_TA, 200);

fprintf('Using window sizes:\n');
fprintf('  a: %d points (%.1f min)\n', window_a, window_a/60);
fprintf('  e: %d points (%.1f min)\n', window_e, window_e/60);
fprintf('  i: %d points (%.1f min)\n', window_i, window_i/60);
fprintf('  Ω: %d points (%.1f min)\n', window_Omega, window_Omega/60);
fprintf('  ω: %d points (%.1f min)\n', window_omega, window_omega/60);
fprintf('  θ: %d points (%.1f min)\n', window_TA, window_TA/60);

% Apply moving average filtering
a_movmean = movmean(keplerian_history(:,1), window_a, "Endpoints", "fill");
e_movmean = movmean(keplerian_history(:,2), window_e, "Endpoints", "fill");
i_movmean = movmean(keplerian_history(:,3), window_i, "Endpoints", "fill");
Omega_movmean = movmean(keplerian_history(:,4), window_Omega, "Endpoints", "fill");
omega_movmean = movmean(keplerian_history(:,5), window_omega, "Endpoints", "fill");
TA_movmean = movmean(keplerian_history(:,6), window_TA, "Endpoints", "fill");

downsample_factor = floor(length(tspan)/10000); % Target ~10,000 points
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

% Now create your plots using the downsampled data if applicable
figure('Name','Geometric Keplerian Elements History - Corrected')
subplot(3,1,1)
hold on
plot(tspan_plot, a_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, a_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Semi-major Axis [km]')
xlabel('Time')
ylabel('a [km]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot, i_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, i_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Inclination [deg]')
xlabel('Time')
ylabel('i [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot, e_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, e_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Eccentricity [-]')
xlabel('Time')
ylabel('e [-]')
legend('Location','best')
grid on

figure('Name','Orientation Keplerian Elements History - Corrected')
subplot(3,1,1)
hold on
plot(tspan_plot, Omega_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, Omega_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Right Ascension Ascending Node [deg]')
xlabel('Time')
ylabel('Ω [deg]')
legend('Location','best')
grid on

subplot(3,1,2)
hold on
plot(tspan_plot, omega_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, omega_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('Argument of Periapsis [deg]')
xlabel('Time')
ylabel('ω [deg]')
legend('Location','best')
grid on

subplot(3,1,3)
hold on
plot(tspan_plot, TA_plot, 'b', 'DisplayName','Propagated', 'LineWidth', 0.5)
plot(tspan_plot, TA_movmean_plot, 'r', 'DisplayName','Filtered', 'LineWidth', 2)
hold off
title('True Anomaly [deg]')
xlabel('Time')
ylabel('θ [deg]')
legend('Location','best')
grid on

%% 3.b. Integration using Gauss Planetary Equations

% Initial Keplerian elements (already in your code)
kep_initial = kep_parameters';  % Convert to column vector: [a; e; i; Omega; omega; theta]

% Create function handle for perturbations
a_per_func = @(t, kep) j2_drag_perturbations_RSW(t, kep, mu_E, J2_E, R_E, omega_E, c_d, AreaOverMass);

% Integration time span (same as before)
tspan_gauss = linspace(0, T_groundtrack, floor(T_groundtrack));

% ODE options
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11, 'MaxStep', 60); % MaxStep = 60s

% Integrate Gauss equations
[TGauss, KepGauss] = ode113(@(t, y) gauss_planetary_equations(t, y, mu_E, a_per_func), ...
                           tspan_gauss, kep_initial, options);

% Convert to Cartesian for comparison
r_gauss = zeros(length(TGauss), 3);
v_gauss = zeros(length(TGauss), 3);

for i = 1:length(TGauss)
    [r_temp, v_temp] = kep2car(KepGauss(i,:)', mu_E);
    r_gauss(i,:) = r_temp';
    v_gauss(i,:) = v_temp';
end

%% 3.c. Compare with Cartesian integration

% Calculate position differences
if length(TPerturbed) == length(TGauss)
    pos_diff = r_perturbed - r_gauss;
    pos_magnitude_diff = sqrt(sum(pos_diff.^2, 2));
    
    figure('Name', 'Comparison: Cartesian vs Gauss Integration');
    subplot(2,1,1);
    plot(TGauss, pos_magnitude_diff);
    xlabel('Time [s]');
    ylabel('Position Difference [km]');
    title('Position Discrepancy between Integration Methods');
    grid on;
    
    subplot(2,1,2);
    semilogy(TGauss, pos_magnitude_diff);
    xlabel('Time [s]');
    ylabel('Position Difference [km] (log scale)');
    title('Logarithmic Scale');
    grid on;
    
    fprintf('Maximum position difference: %.6f km\n', max(pos_magnitude_diff));
    fprintf('Mean position difference: %.6f km\n', mean(pos_magnitude_diff));
end

%% Analyze secular rates from Gauss equations

% Extract filtered secular rates using larger windows
window_secular = floor(points_per_orbit * 10);  % 10 orbits

% Compute rates from Gauss integration
da_dt = gradient(KepGauss(:,1), TGauss);
de_dt = gradient(KepGauss(:,2), TGauss);
di_dt = gradient(KepGauss(:,3), TGauss);
dOmega_dt = gradient(KepGauss(:,4), TGauss);
domega_dt = gradient(KepGauss(:,5), TGauss);

% Filter to get secular rates
da_dt_secular = movmean(da_dt, window_secular, "Endpoints", "fill");
de_dt_secular = movmean(de_dt, window_secular, "Endpoints", "fill");
di_dt_secular = movmean(di_dt, window_secular, "Endpoints", "fill");
dOmega_dt_secular = movmean(dOmega_dt, window_secular, "Endpoints", "fill");
domega_dt_secular = movmean(domega_dt, window_secular, "Endpoints", "fill");

% Plot secular rates
figure('Name', 'Secular Rates from Gauss Equations');
subplot(3,2,1);
plot(TGauss, da_dt_secular * 86400 * 365.25);  % Convert to km/year
xlabel('Time [s]');
ylabel('da/dt [km/year]');
title('Semi-major Axis Rate');
grid on;

subplot(3,2,2);
plot(TGauss, de_dt_secular * 86400 * 365.25);  % Convert to 1/year
xlabel('Time [s]');
ylabel('de/dt [1/year]');
title('Eccentricity Rate');
grid on;

subplot(3,2,3);
plot(TGauss, rad2deg(di_dt_secular) * 86400 * 365.25);  % Convert to deg/year
xlabel('Time [s]');
ylabel('di/dt [deg/year]');
title('Inclination Rate');
grid on;

subplot(3,2,4);
plot(TGauss, rad2deg(dOmega_dt_secular) * 86400 * 365.25);  % Convert to deg/year
xlabel('Time [s]');
ylabel('dΩ/dt [deg/year]');
title('RAAN Rate');
grid on;

subplot(3,2,5);
plot(TGauss, rad2deg(domega_dt_secular) * 86400 * 365.25);  % Convert to deg/year
xlabel('Time [s]');
ylabel('dω/dt [deg/year]');
title('Argument of Perigee Rate');
grid on;

%% Theoretical secular rates from J2 (for comparison)

% Analytical secular rates due to J2 (Vallado, 2013)
n = sqrt(mu_E / a_initial^3);  % Mean motion

% Secular rates (first-order approximation)
dOmega_dt_J2 = -1.5 * J2_E * (R_E^2) * n * cos(i_initial) / (a_initial^2 * (1 - e_initial^2)^2);
domega_dt_J2 = 0.75 * J2_E * (R_E^2) * n * (5*cos(i_initial)^2 - 1) / (a_initial^2 * (1 - e_initial^2)^2);

fprintf('\nTheoretical secular rates from J2:\n');
fprintf('  dΩ/dt: %.6f deg/day\n', rad2deg(dOmega_dt_J2) * 86400);
fprintf('  dω/dt: %.6f deg/day\n', rad2deg(domega_dt_J2) * 86400);

% Compare with numerical results
mean_dOmega_num = mean(dOmega_dt_secular);
mean_domega_num = mean(domega_dt_secular);

fprintf('\nNumerical secular rates (average):\n');
fprintf('  dΩ/dt: %.6f deg/day\n', rad2deg(mean_dOmega_num) * 86400);
fprintf('  dω/dt: %.6f deg/day\n', rad2deg(mean_domega_num) * 86400);

%% 5. Animation

% figure('Name','Animation')
% hold on
% axis equal
% view(3)
% surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none');
% comet3( YPerturbed(:,1), YPerturbed(:,2), YPerturbed(:,3))
% 
% hold off




%% 7. Evolution of two real satellites
% It is over the span of 1 month between 2024-01-01 and 2024-01-28


%TLE Used
%HEO
%0 BREEZE-M DEB (TANK)
%1 37604U 11021C   23365.61578929  .00002286  00000-0  58046-3 0  9995
%2 37604  49.3556 219.8807 5924089 325.1862   7.5827  4.14291330182369

%0 BREEZE-M DEB (TANK)
%1 37604U 11021C   24014.61316955  .00003940  00000-0  10592-2 0  9997
%2 37604  49.3588 211.9309 5920798 332.0673   5.9418  4.14377194182972

%LEO
%0 TITAN 3C TRANSTAGE DEB
%1  1716U 65082BX  23365.87499973  .00000824  00000-0  17701-3 0  9994
%2  1716  32.1887 169.4777 0039702 221.3397 138.4238 14.55334050 77929

%0 TITAN 3C TRANSTAGE DEB
%1 01716U 65082BX  24013.56519027  .00001449  00000-0  31727-3 0  9992
%2 01716  32.1908  95.3036 0039681 333.6628  26.2000 14.55361385 79270

%0 TITAN 3C TRANSTAGE DEB
%1  1716U 65082BX  24028.17581063  .00001062  00000-0  23051-3 0  9996
%2  1716  32.1920   9.9040 0039341 104.2793 256.2224 14.55397974 81907

[utc_time_tank, mission_info_tank] = read_ephemeris_file('Ephemeris_BreezeM_Tank.txt');

[utc_time_titan, mission_info_titan] = read_ephemeris_file('Ephemeris_Titan3C_Transtage.txt');


a_eph_tank = mission_info_tank(:,10);      % Semi-major axis [km]
e_eph_tank = mission_info_tank(:,1);       % Eccentricity (already in column 1)
i_eph_tank = deg2rad(mission_info_tank(:,3));  % Inclination [rad]
OM_eph_tank = deg2rad(mission_info_tank(:,4)); % RAAN [rad]
w_eph_tank = deg2rad(mission_info_tank(:,5));  % Argument of periapsis [rad]
TA_eph_tank = deg2rad(mission_info_tank(:,9)); % True anomaly [rad]

keplerian_elements_eph_tank = [a_eph_tank, e_eph_tank, i_eph_tank, OM_eph_tank, w_eph_tank, TA_eph_tank];
keplerian_elements_eph_tank(:,6) = unwrap(keplerian_elements_eph_tank(:,6));

a_eph_titan = mission_info_titan(:,10);      % Semi-major axis [km]
e_eph_titan = mission_info_titan(:,1);       % Eccentricity (already in column 1)
i_eph_titan = deg2rad(mission_info_titan(:,3));  % Inclination [rad]
OM_eph_titan = deg2rad(mission_info_titan(:,4)); % RAAN [rad]
w_eph_titan = deg2rad(mission_info_titan(:,5));  % Argument of periapsis [rad]
TA_eph_titan = deg2rad(mission_info_titan(:,9)); % True anomaly [rad]

keplerian_elements_eph_titan = [a_eph_titan, e_eph_titan, i_eph_titan, OM_eph_titan, w_eph_titan, TA_eph_titan];

keplerian_elements_eph_titan(:,6) = unwrap(keplerian_elements_eph_titan(:,6));

fprintf('Processing ephemeris ground track for Tank...\n');
gtrack_data_eph_tank = zeros(length(utc_time_tank), 4);

for i = 1:length(utc_time_tank)
    % Convert to Cartesian coordinates
    [r_eph_tank, v_eph_tank] = kep2car(keplerian_elements_eph_tank(i,:), mu_E);
    
    % Time difference from reference
    dt_tank = seconds(utc_time_tank(i) - utc_time_tank(1));
    R_eph_tank = norm(r_eph_tank);
    
    % Calculate declination (latitude in radians)
    delta_eph_tank = asin(r_eph_tank(3) / R_eph_tank);
    gtrack_data_eph_tank(i,1) = delta_eph_tank;
    
    % Calculate right ascension
    alpha_eph_tank = atan2(r_eph_tank(2), r_eph_tank(1));
    gtrack_data_eph_tank(i,2) = alpha_eph_tank;
    
    % Calculate Greenwich sidereal time
    theta_g_eph_tank = theta_g_t0_rad + w_earth_rad * dt_tank;
    
    % Calculate longitude (wrap to [-180, 180])
    lon_eph_tank = rad2deg(alpha_eph_tank - theta_g_eph_tank);
    gtrack_data_eph_tank(i,3) = mod(lon_eph_tank + 180, 360) - 180;
    
    % Calculate latitude in degrees
    gtrack_data_eph_tank(i,4) = rad2deg(delta_eph_tank);
end


% Process ephemeris ground track
fprintf('Processing ephemeris ground track for Titan...\n');
gtrack_data_eph_titan = zeros(length(utc_time_titan), 4);

for i = 1:length(utc_time_titan)
    % Convert to Cartesian coordinates
    [r_eph_titan, v_eph_titan] = kep2car(keplerian_elements_eph_titan(i,:), mu_E);
    
    % Time difference from reference
    dt_titan = seconds(utc_time_titan(i) - utc_time_titan(1));
    R_eph_titan = norm(r_eph_titan);
    
    % Calculate declination (latitude in radians)
    delta_eph_titan = asin(r_eph_titan(3) / R_eph_titan);
    gtrack_data_eph_titan(i,1) = delta_eph_titan;
    
    % Calculate right ascension
    alpha_eph_titan = atan2(r_eph_titan(2), r_eph_titan(1));
    gtrack_data_eph_titan(i,2) = alpha_eph_titan;
    
    % Calculate Greenwich sidereal time
    theta_g_eph_titan = theta_g_t0_rad + w_earth_rad * dt_titan;
    
    % Calculate longitude (wrap to [-180, 180])
    lon_eph_titan = rad2deg(alpha_eph_titan - theta_g_eph_titan);
    gtrack_data_eph_titan(i,3) = mod(lon_eph_titan + 180, 360) - 180;
    
    % Calculate latitude in degrees
    gtrack_data_eph_titan(i,4) = rad2deg(delta_eph_titan);
end

fprintf('Plotting results...\n');

figure('Name','Ground Track for Ephemeris')
imagesc([-180 180], [-90 90], flipud(EarthImage));
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

plot_ground_track(gtrack_data_eph_tank(:,3), gtrack_data_eph_tank(:,4), 'k', '--');
plot_ground_track(gtrack_data_eph_titan(:,3), gtrack_data_eph_titan(:,4), 'b', '.-');

xlabel('Longitude [degrees]')
ylabel('Latitude [degrees]')
title('Satellite Ground Track Comparison')
legend('BREEZE-M DEB (TANK)', 'TITAN 3C TRANSTAGE DEB', 'Location', 'best')
xlim([-180 180])
ylim([-90 90])
grid on
hold off




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