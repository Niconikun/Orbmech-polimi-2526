%% Exercise 1

clc
clear
close all


% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_t1 = [-21800; 37900; 0];
r_t2 = [27300; 27700; 0];

ToF = 15*3600 + 6*60 + 40;

mu_E = astroConstants(13);

[a,P,e,ERROR,v_t1,v_t2,TPAR,THETA] = lambertMR( r_t1, r_t2, ToF, mu_E, 0, 0, 0 );

fprintf('Computed value of a: ');
disp(a)
fprintf('Computed value of v_t1: ');
disp(v_t1)
fprintf('Computed value of v_t2: ');
disp(v_t2)


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


J2 = 0;
R = astroConstants(23);
Rp = 3000;

[X, Y, Z] = sphere(50);
X = R*X;
Y = R*Y;
Z = R*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

[Xp, Yp, Zp] = sphere(50);
Xp = Rp*Xp;
Yp = Rp*Yp;
Zp = Rp*Zp;

S0_fwd = [r_t1; v_t1'];
S0_bwd = [r_t2; v_t2'];

t0 = 0;
tf = 2*pi * sqrt(a^3 / mu_E);
tspan = linspace(t0,tf,20000);




% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Propagation forward in time
[ ~ , S_fwd ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R, J2 ), tspan, S0_fwd, options );

% Propagatino backward in time
[ ~ , S_bwd ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R, J2 ), tspan, S0_bwd, options );




% Plot of the propagation forward in time
figure('Name', 'Exercise 1: propagation forward in time')
grid on
hold on
axis equal
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
plot3(S_fwd(:,1),S_fwd(:,2),S_fwd(:,3),'-r','LineWidth',3)
point1 = surf(Xp + r_t1(1), Yp + r_t1(2), Zp + r_t1(3), 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'P1');
point2 = surf(Xp + r_t2(1), Yp + r_t2(2), Zp + r_t2(3), 'FaceColor', 'g', 'EdgeColor', 'none', 'DisplayName', 'P2');
legend([point1 point2]);
xlabel('x [km]')
ylabel('y [km]')
view(3)


% Plot of the propagation bacward in time
figure('Name', 'Exercise 1: propagation backward in time')
grid on
hold on
axis equal
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
plot3(S_bwd(:,1),S_bwd(:,2),S_bwd(:,3),'-m','LineWidth',3)
point1 = surf(Xp + r_t1(1), Yp + r_t1(2), Zp + r_t1(3), 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'P1');
point2 = surf(Xp + r_t2(1), Yp + r_t2(2), Zp + r_t2(3), 'FaceColor', 'g', 'EdgeColor', 'none', 'DisplayName', 'P2');
legend([point1 point2]);
xlabel('x [km]')
ylabel('y [km]')
view(3)


% As expected the foward propagation and the backward propagation produce
% the same image.


%% Exercise 2

clc
clear
close all

mu_E = astroConstants(13);
R_E  = astroConstants(23);

kep_1 = [12500, 0, 0, 0, 0, 2*pi/3];
kep_2 = [9500, 0.3, 0, 0, 0, 250*pi/180];

ToF = 3300; % [s]


% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r_1, v_1] = kep2car(kep_1, mu_E)
[r_2, v_2] = kep2car(kep_2, mu_E)

S_initial = [r_1; v_1];
S_final = [r_2; v_2];


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,P,e,ERROR,v_t1,v_t2,TPAR,THETA] = lambertMR( r_1, r_2, ToF, mu_E, 0, 0, 0 );

fprintf('Computed value of a: ');
disp(a)
fprintf('Computed value of v_t1: ');
disp(v_t1)
fprintf('Computed value of v_t2: ');
disp(v_t2)


% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta_V1 = norm(v_t1' - v_1)
Delta_V2 = norm(v_2 - v_t2')

Total_Delta_V = Delta_V1 + Delta_V2


% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J2 = 0;
t0 = 0;

T_arc = 2*pi * sqrt(a^3 / mu_E);
tspan_used = linspace(t0,ToF,200000);
tspan_unused = linspace(t0,T_arc-ToF,200000);

T_initial = 2*pi * sqrt(12500^3 / mu_E);
tspan_initial = linspace(t0,T_initial,200000);

T_final = 2*pi * sqrt(9500^3 / mu_E);
tspan_final = linspace(t0,T_final,200000);




S0_used = [r_1; v_t1'];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Propagation of the used transfer arc
[ ~ , S_used ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R_E, J2 ), tspan_used, S0_used, options );

% Propagation of the unused transfer arc
[ ~ , S_unused ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R_E, J2 ), tspan_unused, S_used(end,:), options );

% Propagation of the initial orbit
[ ~ , S_initial ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R_E, J2 ), tspan_initial, S_initial, options );

% Propagation of the initial orbit
[ ~ , S_final ] = ode113( @(t,y) ode_perturbed_2bp( t, y, mu_E, R_E, J2 ), tspan_final, S_final, options );



% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rp = 700;

[X, Y, Z] = sphere(50);
X = R_E*X;
Y = R_E*Y;
Z = R_E*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

[Xp, Yp, Zp] = sphere(50);
Xp = Rp*Xp;
Yp = Rp*Yp;
Zp = Rp*Zp;

% Plot
figure('Name', 'Exercise 2')
grid on
hold on
axis equal

% Earth
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');

% Initial orbit
initial_orbit = plot3( S_initial(:,1), S_initial(:,2), S_initial(:,3), 'Color', [0.2 0.3 0.6], 'LineWidth', 1.5, 'DisplayName','Initial orbit');

% Final orbit
final_orbit = plot3( S_final(:,1), S_final(:,2), S_final(:,3), 'Color', [0.3 0.55 0.3], 'LineWidth', 1.5, 'DisplayName','Final orbit');

% Transfer arc
transfer_arc = plot3( S_used(:,1), S_used(:,2), S_used(:,3), 'Color', [1.0 0.7 0.4], 'LineWidth', 1.5, 'DisplayName', 'Transfer arc');

% Transfer arc
unused_arc = plot3( S_unused(:,1), S_unused(:,2), S_unused(:,3), '--', 'Color', [1.0 0.7 0.4], 'LineWidth', 1.5, 'DisplayName', 'Unused arc');

% Point 1
point1 = surf(Xp + r_1(1), Yp + r_1(2), Zp + r_1(3), 'FaceColor', [0.2 0.3 0.6], 'EdgeColor', 'none', 'DisplayName', 'Transfer initial point');

% Point 2
point2 = surf(Xp + r_2(1), Yp + r_2(2), Zp + r_2(3), 'FaceColor', [0.3 0.55 0.3], 'EdgeColor', 'none', 'DisplayName', 'Transfer final point');

% Other plotting parameters
legend([initial_orbit final_orbit transfer_arc unused_arc point1 point2], 'Location','bestoutside');
xlabel('x [km]')
ylabel('y [km]')
xlim([-2 2] * 1e4)
ylim([-2 2] * 1e4)




%% Exercise 3

clc
clear
close all


% Part 2 and 3 and 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Planetary constant of the Sun
mu_Sun = astroConstants(4);

% Planets involved
Earth = 3;
Mars = 4;

% Initialisation of the window
Departure_time_earliest = [2003, 04, 01, 0, 0, 0.00];
Departure_time_latest   = [2003, 08, 01, 0, 0, 0.00];
Arrival_time_earliest   = [2003, 09, 01, 0, 0, 0.00];
Arrival_time_latest     = [2004, 03, 01, 0, 0, 0.00];

% Conversion to mjd
Departure_time_earliest = date2mjd2000(Departure_time_earliest);
Departure_time_latest   = date2mjd2000(Departure_time_latest);
Arrival_time_earliest   = date2mjd2000(Arrival_time_earliest);
Arrival_time_latest     = date2mjd2000(Arrival_time_latest);



% Resolutions of the departure and arrival windows
res_t1 = 1000;
res_t2 = 1000;

t1 = linspace(Departure_time_earliest, Departure_time_latest, res_t1);
t2 = linspace(Arrival_time_earliest, Arrival_time_latest, res_t2);

Delta_V_tot = zeros(res_t1,res_t2);
ToF = zeros(res_t1,res_t2);

parfor i = 1:res_t1

    for j = 1:res_t2

        Delta_V_tot(j,i) = Lambert_direct_transfer(Earth,t1(i),Mars,t2(j),mu_Sun);
        ToF(j,i) = t2(j) - t1(i);

    end

    fprintf('%.2f %%\n', (1-i/res_t1+1/res_t1)*100)

end



% Refinement of the solution
[DV_col, j] = min(Delta_V_tot);
[DV_min, i] = min(DV_col);
j = j(i);

options = optimoptions('fminunc', 'Display', 'none', 'StepTolerance', 1e-16);

f = @(x) Lambert_transfer_for_fminunc(Earth,x(1),Mars,x(2),mu_Sun);
x0 = [t1(j); t2(i)];
[xmin, fval] = fminunc(f, x0, options);

fprintf('Minimum DV: ');
disp(fval)


% Plot
figure('Name','Exercise 3: 3D plot of the cost')
[X, Y] = meshgrid(t1,t2);
hold on
surf(X,Y,Delta_V_tot,'EdgeColor', 'none')
plot3(xmin(1),xmin(2),fval,'ok','MarkerFaceColor','r','MarkerSize', 5)
view(3)

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);


figure('Name','Exercise 3: Porkchop plot')
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
clabel(C, h);
grid on
hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);



figure('Name','Exercise 3: Porkchop plot with constant ToF')
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
grid on

cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
 
[C2, h2] = contour(X,Y,ToF,[60, 120, 180, 240, 300],'LineColor', 'k');
clabel(C2, h2);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);

% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kep_1, ~] = uplanet (xmin(1), Earth);
[kep_2, ~] = uplanet (xmin(2), Mars);

[r1,v1] = kep2car(kep_1, mu_Sun);
[r2,v2] = kep2car(kep_2, mu_Sun);

Delta_t = (xmin(2) - xmin(1)) * 24 * 3600;

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, Delta_t, mu_Sun, 0, 0, 0 );
kep_t = car2kep(r1,v_t1',mu_Sun);

T_Earth = 365.26 * 24 * 3600;
T_Mars = 2*pi * sqrt( kep_2(1)^3 / mu_Sun);
T_transf = Delta_t;
T_unused = 2*pi *sqrt( kep_t(1)^3 / mu_Sun) - T_transf;

t0 = 0;

tspan_Earth = linspace(t0,T_Earth,20000);
tspan_Mars  = linspace(t0,T_Mars,20000);
tspan_transf = linspace(t0,T_transf,20000);
tspan_unused = linspace(t0,T_unused,20000);

S_Earth = [r1; v1];
S_Mars  = [r2; v2];
S_transf = [r1; v_t1'];
S_unused = [r2; v_t2'];


% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Propagation of Earth orbit
[ ~ , S_Earth ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_Earth, S_Earth, options );

% Propagation of Earth orbit
[ ~ , S_Mars ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_Mars, S_Mars, options );

% Propagation of Earth orbit
[ ~ , S_transf ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_transf, S_transf, options );

% Propagation of Earth orbit
[ ~ , S_unused ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_unused, S_unused, options );


R_T = 15000000;
[X_T, Y_T, Z_T] = sphere(50);
X_T = R_T * X_T;
Y_T = R_T * Y_T;
Z_T = R_T * Z_T;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

R_M = 15000000;
[X_M, Y_M, Z_M] = sphere(50);
X_M = R_M * X_M;
Y_M = R_M * Y_M;
Z_M = R_M * Z_M;
% mars = imread('MarsTexture.jpg');
% mars = flipud(mars);

R_S = 12000000;
[X_S, Y_S, Z_S] = sphere(50);
X_S = R_S * X_S;
Y_S = R_S * Y_S;
Z_S = R_S * Z_S;
% sun = imread('SunTexture.jpg');
% sun = flipud(sun);


figure('Name','Exercise 3: plot of the orbits')

grid on
hold on
axis equal

surf(X_T + r1(1), Y_T + r1(2), Z_T + r1(3), 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
surf(X_M + r2(1), Y_M + r2(2), Z_M + r2(3), 'FaceColor', 'r', 'EdgeColor', 'none');
surf(X_S        , Y_S        , Z_S        , 'FaceColor', 'y', 'EdgeColor', 'y');

plot3(S_Earth(:,1),S_Earth(:,2),S_Earth(:,3),'-b','LineWidth',3)
plot3(S_Mars(:,1) ,S_Mars(:,2), S_Mars(:,3),'-r','LineWidth',3)
plot3(S_transf(:,1),S_transf(:,2),S_transf(:,3),'-g','LineWidth',3)

xlabel('x [km]')
ylabel('y [km]')
view(3)


%% Exercise 4

clc
clear
close all


% Mercury %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

planet = 'Mercury';
Mercury = 1;
Earth = 3;
Departure_time_earliest = [2023 11 01 00 00 00.00];
Departure_time_latest   = [2025 01 01 00 00 00.00];
Arrival_time_earliest   = [2024 04 01 00 00 00.00];
Arrival_time_latest     = [2025 03 01 00 00 00.00];

type = 0;
res_t1 = 200;
res_t2 = 200;

[xmin, fval] = find_best_transfer(Earth, Mercury, Departure_time_earliest, Departure_time_latest, Arrival_time_earliest, Arrival_time_latest, type, res_t1, res_t2, planet);


% Venus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

planet = 'Venus';
Mercury = 2;
Earth = 3;
Departure_time_earliest = [2024 06 01 00 00 00.00];
Departure_time_latest   = [2026 11 01 00 00 00.00];
Arrival_time_earliest   = [2024 12 01 00 00 00.00];
Arrival_time_latest     = [2027 06 01 00 00 00.00];

type = 1;
res_t1 = 200;
res_t2 = 200;

[xmin, fval] = find_best_transfer(Earth, Mercury, Departure_time_earliest, Departure_time_latest, Arrival_time_earliest, Arrival_time_latest, type, res_t1, res_t2, planet);


%% Additional exercise

clc
clear
close all


% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S1_t1 = [-1964.809; 2821.834; 596.808; -2.902; -2.044; 0.107];
S2_t2 = [4836.089; -6945.559; -1468.959; 1.629; 1.147; -0.061];

ToF = 1*3600 + 58.76*60;

mu_M = astroConstants(14);

[a,P,e,ERROR,v_t1,v_t2,TPAR,THETA] = lambertMR( S1_t1(1:3), S2_t2(1:3), ToF, mu_M, 0, 0, 0 );

DV_total = norm(v_t1' - S1_t1(4:6)) + norm(S2_t2(4:6) - v_t2');

fprintf('Computed value of v_t1: ');
disp(v_t1)
fprintf('Computed value of v_t2: ');
disp(v_t2)
fprintf('Computed value of the total cost: ');
disp(DV_total)


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computed value of a: ');
disp(a)


% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kep_initial = car2kep(S1_t1(1:3),S1_t1(4:6),mu_M);
kep_final   = car2kep(S2_t2(1:3),S2_t2(4:6),mu_M);
kep_transf  = car2kep(S1_t1(1:3),v_t1', mu_M);

figure
hold on
grid on
axis equal
plotOrbit (kep_initial, mu_M, 300, 'r')
plotOrbit (kep_final, mu_M, 300, 'g')
plotOrbit (kep_transf, mu_M, 300, 'b')
view(3)
opts.Units = 'km';
planet3D('Mars',opts)

% Sembra un trasferimento bitangente o un Homhann


% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_S = astroConstants(4);

name = 'Mercury -> Venus Transfer';
Mercury = 1;
Venus = 2;
Departure_time_earliest = [2024 06 01 00 00 00.00];
Departure_time_latest   = [2026 11 01 00 00 00.00];
Arrival_time_earliest   = [2024 12 01 00 00 00.00];
Arrival_time_latest     = [2027 06 01 00 00 00.00];

type = 2;
res_t1 = 200;
res_t2 = 200;

[xmin, fval] = find_best_transfer_mod(Mercury, Venus, Departure_time_earliest, Departure_time_latest, Arrival_time_earliest, Arrival_time_latest, type, res_t1, res_t2, name);


[DV_tot, DV_1, DV_2] = Lambert_direct_transfer(1,xmin(1),2,xmin(2),mu_S);


%% FUNCTIONS

function dy = ode_perturbed_2bp( ~, y, mu, R, J2 )
%ode_perturbed_2bp ODE system for the two-body problem (Keplerian motion)
%considering also the perturbation of Earth's gravity field
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% R[1] Radius [L]
% J2[1] second zonal harmonic [-]

%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez, Federico Masiero
%
% VERSIONS
% 2025-10-06: First version
%
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

% Compute the perturbation as an acceleration

a_J2_x = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(1) / rnorm * (5 * r(3)^2 / rnorm^2 - 1));
a_J2_y = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(2) / rnorm * (5 * r(3)^2 / rnorm^2 - 1));
a_J2_z = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(3) / rnorm * (5 * r(3)^2 / rnorm^2 - 3));

a_J2 = [a_J2_x; a_J2_y; a_J2_z];

% Set the derivatives of the state
dy = [ v; 
    (-mu/rnorm^3)*r + a_J2];
end



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




function [DV_tot, DV_1, DV_2] = Lambert_direct_transfer(body1,t1,body2,t2,mu)

ToF = (t2 - t1) * 24 * 3600;

[kep_1, ~] = uplanet (t1, body1);
[kep_2, ~] = uplanet (t2, body2);

[r1,v1] = kep2car(kep_1, mu);
[r2,v2] = kep2car(kep_2, mu);

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, ToF, mu, 0, 0, 0 );

DV_1 = v_t1' - v1;
DV_2 = v2 - v_t2';

DV_tot = norm(DV_1) + norm(DV_2);

end


function DV_tot = Lambert_transfer_for_fminunc(body1,t1,body2,t2,mu)

ToF = (t2 - t1) * 24 * 3600;

[kep_1, ~] = uplanet (t1, body1);
[kep_2, ~] = uplanet (t2, body2);

[r1,v1] = kep2car(kep_1, mu);
[r2,v2] = kep2car(kep_2, mu);

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, ToF, mu, 0, 0, 0 );

DV_1 = v_t1' - v1;
DV_2 = v2 - v_t2';

DV_tot = norm(DV_1) + norm(DV_2);

end




function [xmin, fval] = find_best_transfer(planet1, planet2, Departure_time_earliest, Departure_time_latest, Arrival_time_earliest, Arrival_time_latest, type, res_t1, res_t2, name)

% Planetary constant of the Sun
mu_Sun = astroConstants(4);

% Planets involved
Earth = planet1;
planet = planet2;

% Conversion to mjd
Departure_time_earliest = date2mjd2000(Departure_time_earliest);
Departure_time_latest   = date2mjd2000(Departure_time_latest);
Arrival_time_earliest   = date2mjd2000(Arrival_time_earliest);
Arrival_time_latest     = date2mjd2000(Arrival_time_latest);

t1 = linspace(Departure_time_earliest, Departure_time_latest, res_t1);
t2 = linspace(Arrival_time_earliest, Arrival_time_latest, res_t2);

Delta_V_tot = zeros(res_t1,res_t2);
ToF = zeros(res_t1,res_t2);

parfor i = 1:res_t1

    for j = 1:res_t2

        Delta_V_tot(j,i) = Lambert_direct_transfer(Earth,t1(i),planet,t2(j),mu_Sun);
        ToF(j,i) = t2(j) - t1(i);

    end

    fprintf('%.2f %%\n', (1-i/res_t1+1/res_t1)*100)

end



% Refinement of the solution
[DV_col, j] = min(Delta_V_tot);
[fval, i] = min(DV_col);
j = j(i);

xmin = [t1(i); t2(j)];


if type == 1

    options = optimoptions('fminunc', 'Display', 'none', 'StepTolerance', 1e-16);

    f = @(x) Lambert_transfer_for_fminunc(Earth,x(1),planet,x(2),mu_Sun);

    [xmin, fval] = fminunc(f, xmin, options);

elseif type == 2

    options = optimoptions('fmincon', 'Display', 'none', 'StepTolerance', 1e-16);

    f = @(x) Lambert_transfer_for_fminunc(Earth,x(1),planet,x(2),mu_Sun);

    lb = [xmin(1)-0.01; xmin(2)-0.01];
    ub = [xmin(1)+0.01; xmin(2)+0.01];

    [xmin, fval] = fmincon(f, xmin, [], [], [], [], lb, ub, [], options);           

end


fprintf('Minimum DV: ');
disp(fval)


% Plot
figure('Name',name)
[X, Y] = meshgrid(t1,t2);
hold on
surf(X,Y,Delta_V_tot,'EdgeColor', 'none')
plot3(xmin(1),xmin(2),fval,'ok','MarkerFaceColor','r','MarkerSize', 5)
view(3)

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);


figure('Name',name)
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
clabel(C, h);
grid on
hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);



figure('Name',name)
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
grid on

cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
 
[C2, h2] = contour(X,Y,ToF,[60, 120, 180, 240, 300],'LineColor', 'k');
clabel(C2, h2);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);

% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kep_1, ~] = uplanet (xmin(1), Earth);
[kep_2, ~] = uplanet (xmin(2), planet);

[r1,v1] = kep2car(kep_1, mu_Sun);
[r2,v2] = kep2car(kep_2, mu_Sun);

Delta_t = (xmin(2) - xmin(1)) * 24 * 3600;

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, Delta_t, mu_Sun, 0, 0, 0 );
kep_t = car2kep(r1,v_t1',mu_Sun);

T_Earth = 365.26 * 24 * 3600;
T_planet = 2*pi * sqrt( kep_2(1)^3 / mu_Sun);
T_transf = Delta_t;
T_unused = 2*pi *sqrt( kep_t(1)^3 / mu_Sun) - T_transf;

t0 = 0;

tspan_Earth = linspace(t0,T_Earth,20000);
tspan_planet  = linspace(t0,T_planet,20000);
tspan_transf = linspace(t0,T_transf,20000);
tspan_unused = linspace(t0,T_unused,20000);

S_Earth = [r1; v1];
S_planet  = [r2; v2];
S_transf = [r1; v_t1'];
S_unused = [r2; v_t2'];


% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Propagation of Earth orbit
[ ~ , S_Earth ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_Earth, S_Earth, options );

% Propagation of Earth orbit
[ ~ , S_planet ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_planet, S_planet, options );

% Propagation of Earth orbit
[ ~ , S_transf ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_transf, S_transf, options );

% Propagation of Earth orbit
[ ~ , S_unused ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_unused, S_unused, options );


R_T = 15000000;
[X_T, Y_T, Z_T] = sphere(50);
X_T = R_T * X_T;
Y_T = R_T * Y_T;
Z_T = R_T * Z_T;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

R_M = 15000000;
[X_M, Y_M, Z_M] = sphere(50);
X_M = R_M * X_M;
Y_M = R_M * Y_M;
Z_M = R_M * Z_M;
% mars = imread('MarsTexture.jpg');
% mars = flipud(mars);

R_S = 12000000;
[X_S, Y_S, Z_S] = sphere(50);
X_S = R_S * X_S;
Y_S = R_S * Y_S;
Z_S = R_S * Z_S;
% sun = imread('SunTexture.jpg');
% sun = flipud(sun);


figure('Name',name)

grid on
hold on
axis equal

surf(X_T + r1(1), Y_T + r1(2), Z_T + r1(3), 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
surf(X_M + r2(1), Y_M + r2(2), Z_M + r2(3), 'FaceColor', 'r', 'EdgeColor', 'none');
surf(X_S        , Y_S        , Z_S        , 'FaceColor', 'y', 'EdgeColor', 'y');

plot3(S_Earth(:,1),S_Earth(:,2),S_Earth(:,3),'-b','LineWidth',3)
plot3(S_planet(:,1) ,S_planet(:,2), S_planet(:,3),'-r','LineWidth',3)
plot3(S_transf(:,1),S_transf(:,2),S_transf(:,3),'-g','LineWidth',3)

xlabel('x [km]')
ylabel('y [km]')
view(3)


end

function plotOrbit(PK, mu, n_pti, line, th_i, th_f, r)
% plotOrbit - Plots the orbit relative to the specific parameters
% 
% plotOrbit (PK, mu, n_pti)
% plotOrbit (PK, mu, n_pti, line)
% plotOrbit (PK, mu, n_pti, th_i, th_f)
% plotOrbit (PK, mu, n_pti, line, th_i, th_f)
% plotOrbit (PK, mu, n_pti, line, th_i, th_f, r)
% 
% -------------------------------------------------------------------------
% Input arguments:
% PK            [6x1] Orbital parameters
% mu            [1x1] Gravitational constant
% n_pti         [1x1] Number of point to plot
%                     (a small number of points makes the execution faster)
% line          [-]   If inserted, specifies the line stile 
%                     (could be a string, a 3D RGB vector or a struct)
%                     (struct must be defined as: var_name.color and var_name.style)
% th_i          [1x1] If inserted, plots from th_i to th_f 
% th_f          [1x1] If inserted, plots from th_i to th_f
%                     (without th_i and th_f plots the full orbit)
% r             [3x1] Position vector, if inserted plots the orbit as if
%                     the attaractor is in that coordinates
%
% -------------------------------------------------------------------------
% PK definition:
% PK is a row vector ordered as: [a, e, i , OM, w ,th]
% where the constants are: 
%                                   a    Semi-major axis
%                                   e    Eccentricity
%                                   i    Inclination
%                                   OM   Right Ascension of Ascending Node (RAAN)
%                                   w    Argument of Perifocus
%                                   th   True anomaly
%
% -------------------------------------------------------------------------
% Group 20 [Magarelli, Masiero, Francipane] - v.6.4


% LINE DETERMINATION
if PK(2) > 1
toll = 1e-6;
D = acos(-1/PK(2));
end

flag = 0; % Initial definition

if ischar(line)
    flag = 1;
elseif isstruct(line)
    flag = 2;
elseif isvector(line)
    if (size(line,1)*size(line,2) == 3) && (size(line,1)+size(line,2) == 4)
        flag = 3;
    else
        error('"line" is not a 3D vector')
    end
end

% CODE --------------------------------------------------------------------
switch nargin
    % 3 ARGUMENTS
    case 3

    if PK(2) > 1
        % Theta vector for plotting, in case of Hyperbolic orbit; the asymptote angle is D
        TH = linspace( -D+toll, D-toll, n_pti);
    else
        % Theta vector for plotting, in case of Circular, Elliptic and Parabolic
        TH = linspace( 0, 2*pi, n_pti ); 
    end
    
    R_t = [];
    for k = 1 : length(TH)
        [R_plot_vec] = par2car_speed(PK(1), PK(2), PK(3), PK(4), PK(5), TH(k), mu);
        R_t = [R_t, R_plot_vec];
    end
    
    % Code with 3 arguments
    plot3(R_t(1,:),R_t(2,:),R_t(3,:)) 



    % 4 ARGUMENTS
    case 4

    if PK(2) > 1
        % Theta vector for plotting, in case of Hyperbolic orbit; the asymptote angle is D
        TH = linspace( -D+toll, D-toll, n_pti);
    else
        % Theta vector for plotting, in case of Circular, Elliptic and Parabolic
        TH = linspace( 0, 2*pi, n_pti ); 
    end
    
    R_t = [];
    for k = 1 : length(TH)
        [R_plot_vec] = kep2car([PK(1), PK(2), PK(3), PK(4), PK(5), TH(k)], mu);
        R_t = [R_t, R_plot_vec];
    end

    % Code with 4 arguments
    switch flag
        case 1
            % Flag character
            plot3(R_t(1,:), R_t(2,:), R_t(3,:), line) 
        case 2
            % Flag struct
            plot3(R_t(1,:), R_t(2,:), R_t(3,:), 'Color', line.color, 'LineStyle', line.style)
        case 3
            % Flag vector
            plot3(R_t(1,:),R_t(2,:),R_t(3,:), 'Color', line)
    end



    % 5 ARGUMENTS
    case 5
    % In this case "th_i" is saved in "line" and "th_f" in "th_i"
        th_f = th_i; 
        th_i = line; 
   
    if PK(2) > 1
        % Theta vector for plotting, in case of Hyperbolical orbit
        if (th_i < th_f) && ...
                (   ( (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>2*pi-D) && (th_f<=2*pi) )||...
                ( (th_i>=0) && (th_i<D) && (th_f>=0) && (th_f<D) )||...
                ( (th_i>-D) && (th_i<=0) && (th_f>-D) && (th_f<=0) )   )
            % th in case the initial anomaly is smaller of the final
            TH = linspace( th_i, th_f, n_pti);
        elseif (th_f < th_i) && (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>=0) && (th_f<D)
            % th in case the initial anomaly is smaller of the final
            TH = linspace( th_i-2*pi, th_f, n_pti);
        else
            % th in case one or more anomaly inputs are not valid
            disp('th_i and/or th_f inserted are not valid, full hyperbolic orbit plotted')
            TH = linspace( -D+toll, D-toll, n_pti);
        end
    else
        % Theta vector for plotting, in case of Circular, Elliptical and
        % Parabolic orbit

        % In case th_f < th_i, adds 2*pi to th_f
        if th_f < th_i
            th_f = th_f + 2*pi;
        end
        % In case th_i < th_f
        TH = linspace( th_i, th_f, n_pti );
    end
    
    R_t = [];
    for k = 1 : length(TH)
        [R_plot_vec] = par2car_speed(PK(1), PK(2), PK(3), PK(4), PK(5), TH(k), mu);
        R_t = [R_t, R_plot_vec];
    end
    
    % Code with 5 arguments
    plot3(R_t(1,:), R_t(2,:), R_t(3,:))



    % 6 ARGUMENTS
    case 6

    if PK(2) > 1
        % Theta vector for plotting, in case of Hyperbolical orbit
        if (th_i < th_f) && ...
                (   ( (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>2*pi-D) && (th_f<=2*pi) )||...
                ( (th_i>=0) && (th_i<D) && (th_f>=0) && (th_f<D) )||...
                ( (th_i>-D) && (th_i<=0) && (th_f>-D) && (th_f<=0) )   )
            % th in case the initial anomaly is smaller of the final
            TH = linspace( th_i, th_f, n_pti);
        elseif (th_f < th_i) && (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>=0) && (th_f<D)
            % th in case the initial anomaly is smaller of the final
            TH = linspace( th_i-2*pi, th_f, n_pti);
        else
            % th in case one or more anomaly inputs are not valid
            disp('th_i and/or th_f inserted are not valid, full hyperbolic orbit plotted')
            TH = linspace( -D+toll, D-toll, n_pti);
        end
    else
        % Theta vector for plotting, in case of Circular, Elliptical and
        % Parabolic orbit

        % In case th_f < th_i, adds 2*pi to th_f
        if th_f < th_i
            th_f = th_f + 2*pi;
        end
        % In case th_i < th_f
        TH = linspace( th_i, th_f, n_pti );
    end
    
    R_t = [];
    for k = 1 : length(TH)
        [R_plot_vec] = par2car_speed(PK(1),PK(2),PK(3),PK(4),PK(5),TH(k),mu);
        R_t = [R_t R_plot_vec];
    end
    
    % Code with 6 arguments
    switch flag
        case 1
            % Flag character
            plot3(R_t(1,:), R_t(2,:), R_t(3,:), line) 
        case 2
            % Flag struct
            plot3(R_t(1,:), R_t(2,:), R_t(3,:), 'Color', line.color, 'LineStyle', line.style)
        case 3
            % Flag vector
            plot3(R_t(1,:),R_t(2,:),R_t(3,:), 'Color', line)
    end




    % 7 ARGUMENTS
    case 7
        
        % Verify that r is 3x1 or 1x3 vector
        [sr,sc] = size(r);
        if ~((sr*sc == 3) && (sr+sc == 4))
            error('r assigned is NOT a 3x1 or 1x3 position vector')
        end

        if PK(2) > 1
            % Theta vector for plotting, in case of Hyperbolical orbit
            if (th_i < th_f) && ...
                    (   ( (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>2*pi-D) && (th_f<=2*pi) )||...
                    ( (th_i>=0) && (th_i<D) && (th_f>=0) && (th_f<D) )||...
                    ( (th_i>-D) && (th_i<=0) && (th_f>-D) && (th_f<=0) )   )
                % th in case the initial anomaly is smaller of the final
                TH = linspace( th_i, th_f, n_pti);
            elseif (th_f < th_i) && (th_i>2*pi-D) && (th_i<=2*pi) && (th_f>=0) && (th_f<D)
                % th in case the initial anomaly is smaller of the final
                TH = linspace( th_i-2*pi, th_f, n_pti);
            else
                % th in case one or more anomaly inputs are not valid
                disp('th_i and/or th_f inserted are not valid, full hyperbolic orbit plotted')
                TH = linspace( -D+toll, D-toll, n_pti);
            end
        else
            % Theta vector for plotting, in case of Circular, Elliptical and
            % Parabolic orbit
    
            % In case th_f < th_i, adds 2*pi to th_f
            if th_f < th_i
                th_f = th_f + 2*pi;
            end
            % In case th_i < th_f
            TH = linspace( th_i, th_f, n_pti );
        end

        R_t = [];
        for k = 1 : length(TH)
            [R_plot_vec] = par2car_speed(PK(1),PK(2),PK(3),PK(4),PK(5),TH(k),mu);
            R_t = [R_t R_plot_vec];
        end

        R_t(1,:) = R_t(1,:) + r(1); % Add to plotOrbit2point
        R_t(2,:) = R_t(2,:) + r(2);
        R_t(3,:) = R_t(3,:) + r(3);

        % Code with 7 arguments
        switch flag
            case 1
                % Flag character
                plot3(R_t(1,:), R_t(2,:), R_t(3,:), line) 
            case 2
                % Flag struct
                plot3(R_t(1,:), R_t(2,:), R_t(3,:), 'Color', line.color, 'LineStyle', line.style)
            case 3
                % Flag vector
                plot3(R_t(1,:),R_t(2,:),R_t(3,:), 'Color', line)
        end
end

end
