%% Lab 03

clc
clear
close all

% Constants
muE = astroConstants(13); % [km^3/s^2] Planetary constant of Earth
muM = astroConstants(14); % [km^3/s^2] Planetary constant of Mars
muS = astroConstants(4); % [km^3/s^2] Planetary constant of Sun
Re = astroConstants(23); % [km] medium radius of Earth

% Plot
[X, Y, Z] = sphere(50);
X = Re*X;
Y = Re*Y;
Z = Re*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth); 

%% Ex 1
% Data
r1 = [-21800; 37900; 0]; % [km]
r2 = [27300; 27700; 0]; % [km]
ToF = 15*3600 + 6*60 + 40; % [s]

[a, ~, ~, ~, v1, ~, ~, ~] = lambertMR(r1, r2, ToF, muE, 0, 0, 0, 0);

% Orbit propagation
figure(1); 
hold on;
T = 2*pi*sqrt(a^3/muE); % Orbital period
S0 = [r1; v1'];
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[t, S] = PlotOrbit(S0, [0 T], muE, 0, 'r', options);

plot3(r1(1), r1(2), r1(3), '.b', 'MarkerSize', 50)
plot3(r2(1), r2(2), r2(3), '.g',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('orbit', 'P1', 'P2')
grid on; 
axis equal;

%% Ex 2
kep1 = [12500; 0; 0; 0; 0; 120*pi/180];
kep2 = [9500; 0.3; 0; 0; 0; 250*pi/180];
ToF = 3300; % [s]

[r1, v1] = kep2car(kep1, muE);
[r2, v2] = kep2car(kep2, muE);

[a_T, p_T, e_T, ERROR, v1_T, v2_T, tpar, theta] = lambertMR(r1, r2, ToF, muE, 0, 0, 0, 0);
v1_T = v1_T';
v2_T = v2_T';
dv1 = norm(v1_T - v1);
dv2 = norm(v2 - v2_T);
dv_tot = dv1 + dv2;

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

figure(2)
T1 = 2*pi*sqrt(kep1(1)^3/muE); 
S01 = [r1; v1];
[t1, S1] = PlotOrbit(S01, [0 T1], muE, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muE); 
S0t = [r1; v1_T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muE, 0, 'b', options);
S0to = [r2; v2_T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muE, 0, '--b', options);

T2 = 2*pi*sqrt(kep2(1)^3/muE); 
S02 = [r2; v2];
[t2, S2] = PlotOrbit(S02, [0 T2], muE, 0, 'r', options);

% Plot
hold on;
plot3(r1(1), r1(2), r1(3), '.b', 'MarkerSize', 50)
plot3(r2(1), r2(2), r2(3), '.g',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('orbit 1', 'orbit 2', 'trasnfer arc', 'transfer orbit', 'P1', 'P2')
grid on; 
axis equal;

%% Ex 3
% Mission definition
%  Departure planet: Earth
%  Target planet: Mars
%  Departure time window: 2003 Apr 1 - 2003 Aug 1
%  Arrival time window: 2003 Sep 1 - 2004 Mar 1

dep = [2003, 4, 1, 0, 0, 0;
       2003, 8, 1, 0, 0, 0];
arr = [2003, 9, 1, 0, 0, 0;
       2004, 3, 1, 0, 0, 0];

day_step = 0.5;

levels = [5 6 7 8 9 10];

planets = [3 4];

[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, -1);
best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight in seconds
ToF = (best_arr - best_dep) * 3600 * 24;

% Keplerian element vector
[kep_dep, ksun_dep] = uplanet(best_dep, 3);
[r_dep, v_dep] = kep2car(kep_dep, muS);
[kep_arr_E, ksun_arr_E] = uplanet(best_arr, 3);
[r_arr_E, v_arr_E] = kep2car(kep_arr_E, muS);
[kep_dep_M, ksun_dep_M] = uplanet(best_dep, 4);
[r_dep_M, v_dep_M] = kep2car(kep_dep_M, muS);
[kep_arr, ksun_arr] = uplanet(best_arr, 4);
[r_arr, v_arr] = kep2car(kep_arr, muS);

[a, p, e, ERROR, v1_dep, v2_arr, tpar, theta] = lambertMR(r_dep, r_arr, ToF, muS, 0, 0, 0, 0);
v1_dep = v1_dep';
v2_arr = v2_arr';

% Propagation
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_dep(1)^3/muS); 
S01 = [r_dep; v_dep];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a^3/muS); 
S0t = [r_dep; v1_dep];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_arr; v2_arr];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_arr(1)^3/muS); 
S02 = [r_arr; v_arr];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

% Plot
hold on;
plot3(r_dep(1), r_dep(2), r_dep(3), '.g', 'MarkerSize', 50)
plot3(r_arr_E(1), r_arr_E(2), r_arr_E(3), '.g',  'MarkerSize', 50)
plot3(r_dep_M(1), r_dep_M(2), r_dep_M(3), '.r',  'MarkerSize', 50)
plot3(r_arr(1), r_arr(2), r_arr(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Earth orbit', 'Mars orbit', 'Transfer arc', 'Transfer orbit', 'Earth at dep', 'Earth at arr', 'Mars at dep', 'Mars at dep', 'Sun')
grid on; 
axis equal;

%% Ex 4.1
% Mission definition
%  Departure planet: Earth
%  Target planet: Venus
%  Departure time window: 2024 Jun 01 - 2026 Nov 01
%  Arrival time window: 2024 Dec 01 - 2027 Jun 01

dep = [2024, 06, 01, 0, 0, 0;
       2026, 11, 01, 0, 0, 0];
arr = [2024, 12, 01, 0, 0, 0;
       2027, 06, 01, 0, 0, 0];

day_step = 0.5;
levels = [5 6 7 8 9 10];
planets = [3; 2];

[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, -1);
best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight
ToF = (best_dep - best_arr) * 24 * 3600; % [s]

% Keplerian element vector
[kep_1, ksun_1] = uplanet(best_dep, 3);
[r_1, v_1] = kep2car(kep_1, muS);
[kep_2_E, ksun_2_E] = uplanet(best_arr, 3);
[r_2_E, v_2_E] = kep2car(kep_2_E, muS);
[kep_1_V, ksun_1_V] = uplanet(best_dep, 4);
[r_1_V, v_1_V] = kep2car(kep_1_V, muS);
[kep_2, ksun_2] = uplanet(best_arr, 4);
[r_2, v_2] = kep2car(kep_2, muS);

[a_T, p_T, e_T, ERROR, v_1T, v_2T, tpar, theta] = lambertMR(r_1, r_2, ToF, muS, 0, 0, 0, 0);
v_1T = v_1T';
v_2T = v_2T';

% Propagation and plot
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_1(1)^3/muS); 
S01 = [r_1; v_1];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muS); 
S0t = [r_1; v_1T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_2; v_2T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_2(1)^3/muS); 
S02 = [r_2; v_2];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

plot3(r_1(1), r_1(2), r_1(3), '.g', 'MarkerSize', 50)
plot3(r_2_E(1), r_2_E(2), r_2_E(3), '.g',  'MarkerSize', 50)
plot3(r_1_V(1), r_1_V(2), r_1_V(3), '.r',  'MarkerSize', 50)
plot3(r_2(1), r_2(2), r_2(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Earth orbit', 'Mars orbit', 'Transfer arc', 'Transfer orbit', 'Earth at dep', 'Earth at arr', 'Mars at dep', 'Mars at arr', 'Sun')
grid on; 
axis equal;

%% Ex 4.2
% Mission definition
%  Departure planet: Earth
%  Target planet: Venus
%  Departure time window: 2024 Jun 01 - 2026 Nov 01
%  Arrival time window: 2024 Dec 01 - 2027 Jun 01
%  v_inf = 3.0 km/s

v_inf = 3; % [km/s]
[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, v_inf);
best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight
ToF = (best_dep - best_arr) * 24 * 3600; % [s]

% Keplerian element vector
[kep_1, ksun_1] = uplanet(best_dep, 3);
[r_1, v_1] = kep2car(kep_1, muS);
[kep_2_E, ksun_2_E] = uplanet(best_arr, 3);
[r_2_E, v_2_E] = kep2car(kep_2_E, muS);
[kep_1_V, ksun_1_V] = uplanet(best_dep, 4);
[r_1_V, v_1_V] = kep2car(kep_1_V, muS);
[kep_2, ksun_2] = uplanet(best_arr, 4);
[r_2, v_2] = kep2car(kep_2, muS);

[a_T, p_T, e_T, ERROR, v_1T, v_2T, tpar, theta] = lambertMR(r_1, r_2, ToF, muS, 0, 0, 0, 0);
v_1T = v_1T';
v_2T = v_2T';

% Propagation and plot
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_1(1)^3/muS); 
S01 = [r_1; v_1];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muS); 
S0t = [r_1; v_1T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_2; v_2T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_2(1)^3/muS); 
S02 = [r_2; v_2];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

plot3(r_1(1), r_1(2), r_1(3), '.g', 'MarkerSize', 50)
plot3(r_2_E(1), r_2_E(2), r_2_E(3), '.g',  'MarkerSize', 50)
plot3(r_1_V(1), r_1_V(2), r_1_V(3), '.r',  'MarkerSize', 50)
plot3(r_2(1), r_2(2), r_2(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Earth orbit', 'Mars orbit', 'Transfer arc', 'Transfer orbit', 'Earth at dep', 'Earth at arr', 'Mars at dep', 'Mars at arr', 'Sun')
grid on; 
axis equal;

%% Ex 4.3
% Mission definition
%  Departure planet: Earth
%  Target planet: Venus
%  Departure time window: 2024 Jun 01 - 2026 Nov 01
%  Arrival time window: 2024 Dec 01 - 2027 Jun 01

dep = [2025, 08, 01, 0, 0, 0;
       2031, 01, 01, 0, 0, 0];
arr = [2026, 01, 01, 0, 0, 0;
       2032, 01, 01, 0, 0, 0];

day_step = 0.5;
levels = [5 6 7 8 9 10];
planets = [3; 4];

[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, -1);
best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight
ToF = (best_dep - best_arr) * 24 * 3600; % [s]

% Keplerian element vector
[kep_1, ksun_1] = uplanet(best_dep, 3);
[r_1, v_1] = kep2car(kep_1, muS);
[kep_2_E, ksun_2_E] = uplanet(best_arr, 3);
[r_2_E, v_2_E] = kep2car(kep_2_E, muS);
[kep_1_V, ksun_1_V] = uplanet(best_dep, 4);
[r_1_V, v_1_V] = kep2car(kep_1_V, muS);
[kep_2, ksun_2] = uplanet(best_arr, 4);
[r_2, v_2] = kep2car(kep_2, muS);

[a_T, p_T, e_T, ERROR, v_1T, v_2T, tpar, theta] = lambertMR(r_1, r_2, ToF, muS, 0, 0, 0, 0);
v_1T = v_1T';
v_2T = v_2T';

% Propagation and plot
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_1(1)^3/muS); 
S01 = [r_1; v_1];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muS); 
S0t = [r_1; v_1T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_2; v_2T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_2(1)^3/muS); 
S02 = [r_2; v_2];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

plot3(r_1(1), r_1(2), r_1(3), '.g', 'MarkerSize', 50)
plot3(r_2_E(1), r_2_E(2), r_2_E(3), '.g',  'MarkerSize', 50)
plot3(r_1_V(1), r_1_V(2), r_1_V(3), '.r',  'MarkerSize', 50)
plot3(r_2(1), r_2(2), r_2(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Earth orbit', 'Mars orbit', 'Transfer arc', 'Transfer orbit', 'Earth at dep', 'Earth at arr', 'Mars at dep', 'Mars at arr', 'Sun')
grid on; 
axis equal;

%% Ex 4.4
% Mission definition
%  Departure planet: Earth
%  Target planet: Venus
%  Departure time window: 2024 Jun 01 - 2026 Nov 01
%  Arrival time window: 2024 Dec 01 - 2027 Jun 01
%  v_inf = 3.0 km/s

v_inf = 3.5; % [km/s]
[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, v_inf);
best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight
ToF = (best_dep - best_arr) * 24 * 3600; % [s]

% Keplerian element vector
[kep_1, ksun_1] = uplanet(best_dep, 3);
[r_1, v_1] = kep2car(kep_1, muS);
[kep_2_E, ksun_2_E] = uplanet(best_arr, 3);
[r_2_E, v_2_E] = kep2car(kep_2_E, muS);
[kep_1_V, ksun_1_V] = uplanet(best_dep, 4);
[r_1_V, v_1_V] = kep2car(kep_1_V, muS);
[kep_2, ksun_2] = uplanet(best_arr, 4);
[r_2, v_2] = kep2car(kep_2, muS);

[a_T, p_T, e_T, ERROR, v_1T, v_2T, tpar, theta] = lambertMR(r_1, r_2, ToF, muS, 0, 0, 0, 0);
v_1T = v_1T';
v_2T = v_2T';

% Propagation and plot
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_1(1)^3/muS); 
S01 = [r_1; v_1];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muS); 
S0t = [r_1; v_1T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_2; v_2T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_2(1)^3/muS); 
S02 = [r_2; v_2];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

plot3(r_1(1), r_1(2), r_1(3), '.g', 'MarkerSize', 50)
plot3(r_2_E(1), r_2_E(2), r_2_E(3), '.g',  'MarkerSize', 50)
plot3(r_1_V(1), r_1_V(2), r_1_V(3), '.r',  'MarkerSize', 50)
plot3(r_2(1), r_2(2), r_2(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Earth orbit', 'Mars orbit', 'Transfer arc', 'Transfer orbit', 'Earth at dep', 'Earth at arr', 'Mars at dep', 'Mars at arr', 'Sun')
grid on; 
axis equal;

%% Additional exercise
% Planet: Mars
r1 = [-1964.809; 2821.834; 596.808];
v1 = [-2.902; -2.044; 0.107];
r2 = [4836.089; -6945.559; -1468.959];
v2 = [1.629; 1.147; -0.061];

ToF = 1*3600 + 58.76*60; % [s]
[a, p, e, ERROR, v1_T, v2_T, t_par, theta] = lambertMR(r1, r2, ToF, muM, 0, 0, 0, 0);
v1_T = v1_T';
v2_T = v2_T';

dv1 = norm(v1_T - v1);
dv2 = norm(v2 - v2_T);

dv_tot = dv1 + dv2;

% 1. dv_tot
% another way: [dv_tot, dv1, dv2] = VelocityCost(r1, v1, r2, v2, ToF, muM);
% (same result)
disp(dv_tot)

% 2. a
disp(a)

% 4. dv_tot constrained
dep_1 = [2024, 06, 01, 0, 0, 0];
dep_2 = [2026, 11, 01, 0, 0, 0];

arr_1 = [2024, 12, 01, 0, 0, 0];
arr_2 = [2027, 06, 01, 0, 0, 0];

dep = [dep_1; dep_2];
arr = [arr_1; arr_2];
day_step = 1;
planets = [1; 2];
levels = 50; 
v_inf = 7;

[DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, v_inf);

best_dep = best_date(1);
best_arr = best_date(2);

% Best solution
disp(minDV)
disp(mjd20002date(best_dep))
disp(mjd20002date(best_arr))

% Time of flight
format long g
ToF = (best_arr - best_dep); % [days]
disp(ToF)
ToF = ToF * 24 * 60 * 60;

% Keplerian element vector
[kep_1, ksun_1] = uplanet(best_dep, 1);
[r_1, v_1] = kep2car(kep_1, muS);
[kep_2_E, ksun_2_E] = uplanet(best_arr, 1);
[r_2_E, v_2_E] = kep2car(kep_2_E, muS);
[kep_1_V, ksun_1_V] = uplanet(best_dep, 2);
[r_1_V, v_1_V] = kep2car(kep_1_V, muS);
[kep_2, ksun_2] = uplanet(best_arr, 2);
[r_2, v_2] = kep2car(kep_2, muS);

[a_T, p_T, e_T, ERROR, v_1T, v_2T, tpar, theta] = lambertMR(r_1, r_2, ToF, muS, 0, 0, 0, 0);
v_1T = v_1T';
v_2T = v_2T';

% Propagation and plot
figure()
hold on
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

T1 = 2*pi*sqrt(kep_1(1)^3/muS); 
S01 = [r_1; v_1];
[t1, S1] = PlotOrbit(S01, [0 T1], muS, 0, 'g', options);

Tt = 2*pi*sqrt(a_T^3/muS); 
S0t = [r_1; v_1T];
[tt, St] = PlotOrbit(S0t, [0 ToF], muS, 0, 'b', options);
S0to = [r_2; v_2T];
[tto, Sto] = PlotOrbit(S0to, [ToF Tt], muS, 0, '--b', options);

T2 = 2*pi*sqrt(kep_2(1)^3/muS); 
S02 = [r_2; v_2];
[t2, S2] = PlotOrbit(S02, [0 T2], muS, 0, 'r', options);

plot3(r_1(1), r_1(2), r_1(3), '.g', 'MarkerSize', 50)
plot3(r_2_E(1), r_2_E(2), r_2_E(3), '.g',  'MarkerSize', 50)
plot3(r_1_V(1), r_1_V(2), r_1_V(3), '.r',  'MarkerSize', 50)
plot3(r_2(1), r_2(2), r_2(3), '.r',  'MarkerSize', 50)
plot3(0, 0, 0, '.y',  'MarkerSize', 50)
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
view(3)
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]');
legend('Mercury orbit', 'Venus orbit', 'Transfer arc', 'Transfer orbit', 'Mercury at dep', 'Mercury at arr', 'Venus at dep', 'Venus at arr', 'Sun')
grid on; 
axis equal;

%% Functions

function dydt = ode_2bp(~, y, muE)

    % position
    r = y(1:3);
    r_norm = norm(r);

    % velocity
    v = y(4:6);

    % acceleration
    a = - (muE / (r_norm)^3) * r;

    dydt = [v; a];

end
