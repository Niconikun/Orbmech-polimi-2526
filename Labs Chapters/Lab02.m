%% Lab02
clc 
clear 
close all 

mu = astroConstants(13); 
J2 = astroConstants(9); 
R = astroConstants(23); 
omega_E = 15.04*pi/180/3600; % [deg/s] 
T_e = 24 * 60 * 60; 
lon_G0 = 0; 

%% Case 1 
% Orbit parameters 
a = 8350; % [km] 
e = 0.1976; % [-] 
i = 60*pi/180; % [rad] 
Omega = 270*pi/180; % [rad] 
omega = 45*pi/180; % [rad] 
f0 = 230*pi/180; % [rad] 
r0 = [-4578.219; -801.084; -7929.708]; % [km] 
v0 = [0.800; -6.0387; 1.385]; % [km/s] 
S0 = [r0; v0]; 

% Orbital elements
h = cross(r0, v0); 
e = cross(v0, h) / mu - r0 / norm(r0); 
a = (norm(h)^2) / ( mu * (1 -norm(e)^2 ) ); 
T = 2 * pi * sqrt(a^3/mu); 

% Time vector
t_vec = linspace(0, 15*T, 10000); 

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_1 = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '.r';
[alpha, delta, lon, lat] = groundTrack(r_1, lon_G0, omega_E, t_vec, str);

% Reapting ground tracks
k = 12;
m = 1;

% Orbit
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0.1976; 
kep = [a; e; i; Omega; omega; f0];
[r0, v0] = kep2car(kep, mu);
S0 = [r0; v0];

% Time vector
T = 2 * pi * sqrt(a^3/mu);
t_vec = linspace(0, 15*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_2 = S(:, 1:3); 

% Ground Track
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r_2, lon_G0, omega_E, t_vec, str);

%% Case 2
% Orbit parameters
a = 26600; % [km]
e = 0.74; % [-]
i = 63.4*pi/180; % [rad]
Omega = 50*pi/180; % [rad]
omega = 280*pi/180; % [rad]
f0 = 0*pi/180; % [rad]
r0 = [3108.128; -1040.299; -6090.022]; % [km]
v0 = [5.743; 8.055; 1.555]; % [km/s]
S0 = [r0; v0];

% Orbital elements
h = cross(r0, v0); 
e = cross(v0, h) / mu - r0 / norm(r0); 
a = (norm(h)^2) / ( mu * (1 -norm(e)^2 ) ); 
T = 2 * pi * sqrt(a^3/mu);

% Time vector
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_3 = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '*r';
[alpha, delta, lon, lat] = groundTrack(r_3, lon_G0, omega_E, t_vec, str);

% Reapting ground tracks
k = 2;
m = 1;

% Orbit
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0.74;
kep = [a; e; i; Omega; omega; f0];
[r0, v0] = kep2car(kep, mu);
S0 = [r0; v0];

% Time vector
T = 2 * pi * sqrt(a^3/mu);
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_4 = S(:, 1:3); 

% Ground Track
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r_4, lon_G0, omega_E, t_vec, str);

%% Case 3
% Orbit parameters
a = 800 + R; % [km]
e = 0; % [-]
i = 0*pi/180; % [rad]
Omega = 0*pi/180; % [rad]
omega = 40*pi/180; % [rad]
f0 = 0*pi/180; % [rad]
r0 = [5493.312; 4609.436; 0.000]; % [km]
v0 = [-4.792; 5.711; 0.000]; % [km/s]
S0 = [r0; v0];

% Orbital elements
h = cross(r0, v0); 
e = cross(v0, h) / mu - r0 / norm(r0); 
a = (norm(h)^2) / ( mu * (1 -norm(e)^2 ) ); 
T = 2 * pi * sqrt(a^3/mu);

% Time vector
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_5 = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '*r';
[alpha, delta, lon, lat] = groundTrack(r_5, lon_G0, omega_E, t_vec, str);

% Reapting ground tracks
k = 12;
m = 1;

% Orbit
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0; 
kep = [a; e; i; Omega; omega; f0];
[r0, v0] = kep2car(kep, mu);
S0 = [r0; v0];

% Time vector
T = 2 * pi * sqrt(a^3/mu);
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_6 = S(:, 1:3); 

% Ground Track
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r_6, lon_G0, omega_E, t_vec, str);

%% Case 4
% Orbit parameters
a = 800 + R; % [km]
e = 0; % [-]
i = 30*pi/180; % [rad]
Omega = 0*pi/180; % [rad]
omega = 40*pi/180; % [rad]
f0 = 0*pi/180; % [rad]
r0 = [5493.312; 3991.889; 2304.718]; % [km]
v0 = [-4.792; 4.946; 2.856]; % [km/s]
S0 = [r0; v0];
lon_G0 = 0;

% Orbital elements
h = cross(r0, v0); 
e = cross(v0, h) / mu - r0 / norm(r0); 
a = (norm(h)^2) / ( mu * (1 -norm(e)^2 ) ); 
T = 2 * pi * sqrt(a^3/mu);

% Time vector
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_7 = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '*r';
[alpha, delta, lon, lat] = groundTrack(r_7, lon_G0, omega_E, t_vec, str);

% Repeating ground tracks
k = 12;
m = 1;

% Orbit
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0; 
kep = [a; e; i; Omega; omega; f0];
[r0, v0] = kep2car(kep, mu);
S0 = [r0; v0];

% Time vector
T = 2 * pi * sqrt(a^3/mu);
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_8 = S(:, 1:3); 

% Ground Track
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r_8, lon_G0, omega_E, t_vec, str);

%% Case 5
% Orbit parameters
a = 800 + R; % [km]
e = 0; % [-]
i = 98*pi/180; % [rad]
Omega = 0*pi/180; % [rad]
omega = 40*pi/180; % [rad]
f0 = 0*pi/180; % [rad]
r0 = [5493.312; -641.510; 4564.578]; % [km]
v0 = [-4.792; -0.795; 5.656]; % [km/s]
S0 = [r0; v0];
lon_G0 = 0;

% Orbital elements
h = cross(r0, v0); 
e = cross(v0, h) / mu - r0 / norm(r0); 
a = (norm(h)^2) / ( mu * (1 -norm(e)^2 ) ); 
T = 2 * pi * sqrt(a^3/mu);

% Time vector
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_9 = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '*r';
[alpha, delta, lon, lat] = groundTrack(r_9, lon_G0, omega_E, t_vec, str);

% Reapting ground tracks
k = 12;
m = 1;

% Orbit
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0; 
kep = [a; e; i; Omega; omega; f0];
[r0, v0] = kep2car(kep, mu);
S0 = [r0; v0];

% Time vector
T = 2 * pi * sqrt(a^3/mu);
t_vec = linspace(0, 30*T, 10000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r_10 = S(:, 1:3); 

% Ground Track
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r_10, lon_G0, omega_E, t_vec, str);

%% Test of car2kep and kep2car
% Case 1 (ok)

% Generic orbit
r0 = [-4578.219; -801.084; -7929.708];  % [km]
v0 = [0.800; -6.0387; 1.385]; % [km/s]
kep = car2kep(r0, v0, mu)

% a = 8350; % [km]
% e = 0.1976; % [-]
% i = 60*pi/180; % [rad] 1.0472
% Omega = 270*pi/180; % [rad] 4.7124
% omega = 45*pi/180; % [rad] 0.7854
% f0 = 230*pi/180; % [rad] 4.0143

%% Case 2 (ok)
% a = 26600; % [km]
% e = 0.74; % [-]
% i = 63.4*pi/180; % [rad] 1.1065
% Omega = 50*pi/180; % [rad] 0.8727
% omega = 280*pi/180; % [rad] 4.8869
% f0 = 0*pi/180; % [rad]

r0 = [3108.128; -1040.299; -6090.022]; % [km]
v0 = [5.743; 8.055; 1.555]; % [km/s]
kep = car2kep(r0, v0, mu)

%% Case 3 (ok)
% a = 800 + R; % [km] 7.1710e+03
% e = 0; % [-]
% i = 0*pi/180; % [rad]
% Omega = 0*pi/180; % [rad]
% omega = 40*pi/180; % [rad]
% f0 = 0*pi/180; % [rad]

r0 = [5493.312; 4609.436; 0.000]; % [km]
v0 = [-4.792; 5.711; 0.000]; % [km/s]
kep = car2kep(r0, v0, mu)

%% Case 4 (ok)
% a = 800 + R; % [km] 7.1710e+03
% e = 0; % [-]
% i = 30*pi/180; % [rad] 0.5236
% Omega = 0*pi/180; % [rad]
% omega = 40*pi/180; % [rad] 0.6981
% f0 = 0*pi/180; % [rad]

r0 = [5493.312; 3991.889; 2304.718]; % [km]
v0 = [-4.792; 4.946; 2.856]; % [km/s]
kep = car2kep(r0, v0, mu)

%% Case 5 (ok)
% a = 800 + R; % [km] 7.1710e+03
% e = 0; % [-]
% i = 98*pi/180; % [rad] 1.7104
% Omega = 0*pi/180; % [rad]
% omega = 40*pi/180; % [rad] 0.6981
% f0 = 0*pi/180; % [rad]

r0 = [5493.312; -641.510; 4564.578]; % [km]
v0 = [-4.792; -0.795; 5.656]; % [km/s]
kep = car2kep(r0, v0, mu)

%% LEMU_NGE
% TLE:
%  1 60532U 24149BS  25317.19011302  .00008386  00000-0  76017-3 0  9990
%  2 60532  97.7064  31.2709 0001859  51.7358 308.4030 14.95586792 67686

% Extraction of data
filename = 'LEMUNGE_car_horizons_results.txt';
[M] = readHorizonsCar(filename);
r = [M(:, 1), M(:, 2), M(:, 3)];

% Time
t_vec = 0 : 2*60 : 2*60*(size(r, 1)-1);

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Longitude of Greenwich meridian
% 2025-Nov-13 00:00:00.0000
date = [2025, 11, 13, 0, 0, 0];
lon_G0 = Greenwich_longitude(date);
lon_G0 = lon_G0 * pi / 180;

% Ground Track 1
str = '.r';
[alpha, delta, lon, lat] = groundTrack(r, lon_G0, omega_E, t_vec, str);

% Estraction of the data
r0 = [M(1, 1); M(1, 2); M(1, 3)];
v0 = [M(1, 4); M(1, 5); M(1, 6)];
S0 = [r0; v0];

% Time
t_vec = 0 : 60 : 60*(size(r, 1)-1);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r = S(:, 1:3); 

% Ground Track 2
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r, lon_G0, omega_E, t_vec, str);

%% ISS
% Extraction of data
filename = 'ISS_results.txt';
[M] = readHorizonsCar(filename);
r = [M(:, 1), M(:, 2), M(:, 3)];

% Time
t_vec = linspace(0, 24*3600, 8641);

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Longitude of Greenwich meridian
date = [2023, 11, 14, 0, 0, 0];
lon_G0 = Greenwich_longitude(date);
lon_G0 = lon_G0 * pi / 180;

% Ground Track 1
str = '.r';
[alpha, delta, lon, lat] = groundTrack(r, lon_G0, omega_E, t_vec, str);

% Estraction of the data
r0 = [M(1, 1); M(1, 2); M(1, 3)];
v0 = [M(1, 4); M(1, 5); M(1, 6)];
S0 = [r0; v0];

% Time
t_vec = linspace(0, 24*3600, 8641);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, 0, R), t_vec, S0, options); 
r = S(:, 1:3); 

% Ground Track 2
str = '.g';
[alpha, delta, lon, lat] = groundTrack(r, lon_G0, omega_E, t_vec, str);

%% ADDITIONAL EXERCISE
% Data
a = 11376; % [km]
e = 0.4178; % [-]
i = 66*pi/180; % [rad] 1.1519
Omega = 40*pi/180; % [rad] 0.6981
omega = 20*pi/180; % [rad] 0.3491
theta0 = 23*pi/180; % [rad] 0.4014
kep0 = [a; e; i; Omega; omega; theta0];
lon_G0 = 0; % [deg]
t0 = 0; % [s]

% 1.1 Cartesian state
[r0, v0] = kep2car(kep0, mu);
disp("1. Cartesian state vector at time 0")
disp(r0)
disp(v0)

% 1.2 Validation of the transformation (ok)
kep_val = car2kep(r0, v0, mu);

% 2.1 Propagation in 2BP
S0 = [r0; v0];
tf = 3 * 24 * 60 * 60;
t_vec = linspace(t0, tf, 1000);

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, 0, R), t_vec, S0, options); 
r = S(:, 1:3); 

% Earth Texture
img = imread('EarthTexture.jpg');
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Ground Track
str = '.r';
[alpha_1, delta_1, lon_1, lat_1] = groundTrack(r, lon_G0, omega_E, t_vec, str);

% 2.2 Propagation in 2BP-perturbed
% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r = S(:, 1:3); 

% Ground Track
str = '*g';
[alpha_2, delta_2, lon_2, lat_2] = groundTrack(r, lon_G0, omega_E, t_vec, str);

% 2.3 Absolute difference between lon and lat
lon_diff = abs(lon_1(end) - lon_2(end))
lat_diff = abs(lat_1(end) - lat_2(end))
disp("2. Absolute difference between latitude and longitude in 2BP perturber and 2BP unperturbed")
disp(lat_diff)
disp(lon_diff)

% 3.1 New orbit for repeating ground track
k = 9;
m = 3;
a = RepeatingGroundTrack(k, m, omega_E, mu);
e = 0.4178; % [-]
i = 66*pi/180; % [rad] 1.1519
Omega = 40*pi/180; % [rad] 0.6981
omega = 20*pi/180; % [rad] 0.3491
theta0 = 23*pi/180; % [rad] 0.4014

kep_3 = [a; e; i; Omega; omega; theta0];
[r0, v0] = kep2car(kep_3, mu);
S0 = [r0; v0];

% Integration
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t_vec, S] = ode45(@(t, y) ode_2bp(t, y, mu, J2, R), t_vec, S0, options); 
r = S(:, 1:3); 

% 3.2 Longitude and latitude after 3 days
str = '.k';
[alpha_3, delta_3, lon_3, lat_3] = groundTrack(r, lon_G0, omega_E, t_vec, str);

disp("3. latitude and longitude after 3 days of orbital evolution in 2BP")
disp(lat_3(end))
disp(lon_3(end))

% 4.1 Cartesian State Vector
% COSMOS 249 DEB (NORAD ID: 3509) 10 March 2025 00:00:00.0000
% TLE: 
%   1  3509U 68091F   25069.54529715  .00083545  00000-0  19788-2 0  9991
%   2  3509  62.3289 188.2731 0159710  65.6486 296.1200 15.29049440836990

filename = '3509_results.txt';
[M] = readHorizonsCar(filename);
r = [M(end, 1); M(end, 2); M(end, 3)];
v = [M(end, 4); M(end, 5); M(end, 6)];

disp("4. Cartesian state vector (10 March 2025 00:00:00.0000) given by Horizon system")
disp(r)
disp(v)

r = [M(:, 1), M(:, 2), M(:, 3)];

% Time
t_vec = linspace(0, 24*3600, 8641);

% Earth Texture
img = imread('EarthTexture.jpg');

% Plot of the track
figure('Name', 'Ground Track')
imshow(img, 'XData', [-180 180], 'YData', [90 -90]);
set(gca, 'YDir', 'normal');
hold on

% Longitude of Greenwich meridian
date = [2025, 03, 09, 0, 0, 0];
lon_G0 = Greenwich_longitude(date);
lon_G0 = lon_G0 * pi / 180;

% Ground Track 1
str = '.r';
[alpha, delta, lon, lat] = groundTrack(r, lon_G0, omega_E, t_vec, str);

%% FUNCTION 

function dydt = ode_2bp(~, y, mu, J2, R)

    % Extract position and velocity
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    x = r(1); 
    y = r(2); 
    z = r(3);

    % Gravitational acceleration
    a_grav = -(mu / r_norm^3) * r;

    % J2 perturbation
    factor = 1.5 * J2 * mu * (R^2) / (r_norm^5);
    zx = (5 * z^2 / r_norm^2);
    aJ2x = factor * x * (zx - 1);
    aJ2y = factor * y * (zx - 1);
    aJ2z = factor * z * (zx - 3);
    aJ2 = [aJ2x; aJ2y; aJ2z];

    % Total acceleration
    a_total = a_grav + aJ2;

    % Derivative
    dydt = [v; a_total];

end