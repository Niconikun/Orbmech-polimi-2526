clc
clear all
close all

addpath("ephemerides\");
addpath("time\");
addpath("Functions 2bp\");

AU = 149597870.700; %[km]

mu_Sun = astroConstants(4);

mjd2000 = date2mjd2000([2030,1,1,0,0,0]);

%% Asteroid
[kep, mass, M] = ephAsteroids( mjd2000, 257323);

[r1, v1] = kep2car(kep,mu_Sun);

y_1= [r1; v1];
T_1 = 2*pi*sqrt(kep(1)^3/mu_Sun); % Orbital period [s]
tspan_1 = linspace( 0, T_1, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_1, y_1, options );

% Calculate ground track
r_1 = Y1(:,1:3);
%gtrack_data = zeros(length(T),4);




%% Jupiter

[kep_jupiter, mass] = uplanet( mjd2000, 5);

[r_jupiter, v_jupiter] = kep2car(kep_jupiter,mu_Sun);

y_jupiter= [r_jupiter; v_jupiter];
T_jupiter = 2*pi*sqrt(kep_jupiter(1)^3/mu_Sun); % Orbital period [s]
tspan_jupiter = linspace( 0, T_jupiter, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T_jupiter, Y_jupiter ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_jupiter, y_jupiter, options );

% Calculate ground track
r_jupiter = Y_jupiter(:,1:3);

%% Mars

[kep_mars, mass] = uplanet( mjd2000, 4);

[r_mars, v_mars] = kep2car(kep_mars,mu_Sun);

y_mars= [r_mars; v_mars];
T_mars = 2*pi*sqrt(kep_mars(1)^3/mu_Sun); % Orbital period [s]
tspan_mars = linspace( 0, T_mars, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T_mars, Y_mars ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_mars, y_mars, options );

% Calculate ground track
r_mars = Y_mars(:,1:3);

%% Earth

[kep_earth, mass] = uplanet( mjd2000, 3);

[r_earth, v_earth] = kep2car(kep_earth,mu_Sun);

y_earth= [r_earth; v_earth];
T_earth = 2*pi*sqrt(kep_earth(1)^3/mu_Sun); % Orbital period [s]
tspan_earth = linspace( 0, T_earth, 10000); % Reduced to 10 periods for clarity

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );

% Perform the integration
[ T_earth, Y_earth ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_earth, y_earth, options );

% Calculate ground track
r_earth = Y_earth(:,1:3);


%% Plot
figure(1)

hold on
plot3(r_1(:,1)./AU, r_1(:,2)./AU, r_1(:,3)./AU,'DisplayName', 'Asteroid')
plot3(r_mars(:,1)./AU, r_mars(:,2)./AU, r_mars(:,3)./AU, 'DisplayName', 'Mars')
plot3(r_jupiter(:,1)./AU, r_jupiter(:,2)./AU, r_jupiter(:,3)./AU, 'DisplayName', 'Jupiter')
plot3(r_earth(:,1)./AU, r_earth(:,2)./AU, r_earth(:,3)./AU, 'DisplayName', 'Earth')
grid on;


plot3(r_1(1,1)./AU, r_1(1,2)./AU, r_1(1,3)./AU, 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
plot3(r_jupiter(1,1)./AU, r_jupiter(1,2)./AU, r_jupiter(1,3)./AU, 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
plot3(r_mars(1,1)./AU, r_mars(1,2)./AU, r_mars(1,3)./AU, 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');
plot3(r_earth(1,1)./AU, r_earth(1,2)./AU, r_earth(1,3)./AU, 'ro', 'MarkerSize', 10, 'DisplayName', 'Initial Position');

legend('Location', 'best');
view(45, 30);
xlabel('X Position (AU)');
ylabel('Y Position (AU)');
zlabel('Z Position (AU)');
