clc
clearvars
close all

mu_E = astroConstants(13);      % Sun's gravitational parameter [km^3/s^2];
R_E = astroConstants(23);  % Earth radius [km]
J2_E = astroConstants(33); % Earth J2
textureImage = imread('EarthTexture.jpg');
ToF = 15*3600 + 5*60 + 40;      % Time in [s];
r1 = [-21800; 37900;0];       % Initial position vector [km]
r2 = [27300;27700;0];     % Final position vector [km]

orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 2;

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, mu_E, orbitType, Nrev, Ncase, optionsLMR);
fprintf("Applied Battin's Algorithm...");


%% Get Transfer Orbit

y_transfer = [r1; VI'];
T_period_transfer = 2*pi*sqrt(A^3/mu_E); % Orbital period [s]
tspan_transfer = linspace( 0, T_period_transfer, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ Ttransfer, Ytransfer ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan_transfer, y_transfer, options );

% Calculate ground track
r_transfer = Ytransfer(:,1:3);
%gtrack_data = zeros(length(T),4);

% In your main script, before the ground track loop:
fprintf('First position vector: [%f, %f, %f]\n', r_transfer(1,1), r_transfer(1,2), r_transfer(1,3));
fprintf('Time range: %f to %f seconds\n', tspan_transfer(1), tspan_transfer(end));


fprintf('Transfer orbit calculated.\n');


%% Plot Initial, Final & Transfer 

% Create the sphere
[XS,YS,ZS] = sphere(50);

scaleFactor = 6371; % Radius of Earth in kilometers; adjust as needed

% Scale the sphere
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;

% Plot the results
fprintf("Plotting the results...");


figure(1)

hold on
%plot3(r_initial(:,1), r_initial(:,2), r_initial(:,3))
%plot3(r_final(:,1), r_final(:,2), r_final(:,3))
plot3(r_transfer(:,1),r_transfer(:,2),r_transfer(:,3))
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
legend('Initial Orbit', 'Final Orbit', 'Transfer Orbit');
grid on;
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none');

plot3(r1(1), r1(2), r1(3), 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
plot3(r2(1), r2(2), r2(3), 'go', 'MarkerSize', 10, 'DisplayName', 'Final Position');

axis equal