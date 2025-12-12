clc
clearvars
close all

addpath("ephemerides\")
addpath("Functions 2bp\")
addpath("time\")


mu_E = astroConstants(13);      % Sun's gravitational parameter [km^3/s^2];
R_E = astroConstants(23);  % Earth radius [km]
J2_E = astroConstants(33); % Earth J2
textureImage = imread('EarthTexture.jpg');
ToF = 3300;      % Time in [s];
kep_1 = [12500,0,0,0,0,120];       % Initial position vector [km]
kep_2 = [9500,0.3,0,0,0,250];     % Final position vector [km]


orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 2;

%% Initial position & Orbit Propagation
[r1, v1] = kep2car(kep_1,mu_E);

y_1= [r1; v1];
T_1 = 2*pi*sqrt(kep_1(1)^3/mu_E); % Orbital period [s]
tspan_1 = linspace( 0, T_1, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T1, Y1 ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan_1, y_1, options );

% Calculate ground track
r_1 = Y1(:,1:3);
%gtrack_data = zeros(length(T),4);

% In your main script, before the ground track loop:
fprintf('First position vector: [%f, %f, %f]\n', r_1(1,1), r_1(1,2), r_1(1,3));
fprintf('Time range: %f to %f seconds\n', tspan_1(1), tspan_1(end));



%% Final Position & Orbit Propagation

[r2,v2] = kep2car(kep_2, mu_E);

y_2 = [r2; v2];
T_2 = 2*pi*sqrt(kep_2(1)^3/mu_E); % Orbital period [s]
tspan_2 = linspace( 0, T_2, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T2, Y2 ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan_2, y_2, options );

% Calculate ground track
r_2 = Y2(:,1:3);
%gtrack_data = zeros(length(T),4);

% In your main script, before the ground track loop:
fprintf('First position vector: [%f, %f, %f]\n', r_2(1,1), r_2(1,2), r_2(1,3));
fprintf('Time range: %f to %f seconds\n', tspan_2(1), tspan_2(end));




%% Lambert Solution
[A,P,E,~,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, mu_E, orbitType, Nrev, Ncase, optionsLMR);
fprintf("Applied Battin's Algorithm...");
fprintf('Initial velocity vector: [%f, %f, %f]\n', VI(1), VI(2), VI(3));
fprintf('Final velocity vector: [%f, %f, %f]\n', VF(1), VF(2), VF(3));

y_transfer = [r1; VI'];
T_transfer = 2*pi*sqrt(A.^3/mu_E); % Orbital period [s]

tspan_transfer = linspace(0, T_transfer, 1000); % Reduced to 10 periods for clarity
tspan_tof = linspace(0,ToF,1000);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ Ttransfer, Y_transfer ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan_transfer, y_transfer, options );
[Ttof, Y_tof] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan_tof, y_transfer, options );

% Calculate ground track
r_transfer = Y_transfer(:,1:3);
r_tof = Y_tof(:,1:3);
%gtrack_data = zeros(length(T),4);

% In your main script, before the ground track loop:
fprintf('First position vector: [%f, %f, %f]\n', r_transfer(1,1), r_transfer(1,2), r_transfer(1,3));
fprintf('Time range: %f to %f seconds\n', tspan_transfer(1), tspan_transfer(end));




%% Maneuver Cost

%calculate the cost of the maneuver
cost_1 = norm(VI-v1');
cost_2 = norm(VF-v2');

cost_total = cost_1 + cost_2;
% Display the results of the maneuver
fprintf('1st Maneuver cost: %.3f km/s\n', cost_1);
fprintf('2nd Maneuver cost: %.3f km/s\n', cost_2);
fprintf('Total Maneuver cost: %.3f km/s\n', cost_total);
fprintf('Transfer Angle: %.3f degrees\n', rad2deg(THETA));
fprintf('Semi-Major Axis: %.3f degrees \n', A);


%% Propagation Transfer Arc


%% Plot of Orbits

% Create the sphere
[XS,YS,ZS] = sphere(50);

scaleFactor = R_E; % Radius of Earth in kilometers; adjust as needed

% Scale the sphere
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;

% Plot the results
fprintf("Plotting the results...");


figure(1)

hold on
plot3(r_1(:,1), r_1(:,2), r_1(:,3))
plot3(r_2(:,1), r_2(:,2), r_2(:,3))
plot3(r_transfer(:,1),r_transfer(:,2),r_transfer(:,3))
plot3(r_tof(:,1),r_tof(:,2),r_tof(:,3))
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
legend('Initial Orbit', 'Final Orbit','Unused Transfer Orbit', 'Transfer Orbit');
grid on;
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none', 'DisplayName','Earth');

plot3(r_1(1,1), r_1(1,2), r_1(1,3), 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
plot3(r_2(1,1), r_2(1,2), r_2(1,3), 'go', 'MarkerSize', 10, 'DisplayName', 'Final Position');

axis equal

figure(2)
plot3(r_transfer(:,1),r_transfer(:,2),r_transfer(:,3))