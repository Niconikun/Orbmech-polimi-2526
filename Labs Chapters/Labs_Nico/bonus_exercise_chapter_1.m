%Bonus Exercises

clc
clearvars
addpath("Functions 2bp\");

r0 = [242.9; 5737.8; 5809.1]; %km
v0 = [-5.2090; -3.2007; 3.3786]; % km/s
y0 = [ r0; v0];

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(29);
J2_E = astroConstants(9);



% Question 1

% Set time span
t_sim = 365*24*60*60; 
tspan = linspace( 0, t_sim, 5000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan, y0, options );

r = Y(:,1:3);
v = Y(:,4:6);

h_v = zeros(length(T), 3);
h_module = zeros(length(T), 1);

for i = 1:length(T)
    [h_v(i,:), h_module(i), ~] = angular_momentum(r(i,:), v(i,:));
end


figure(1)
% get specificenergy from all variables of r and v
specificEnergy = zeros(length(T), 1);
for i = 1:length(T)
    specificEnergy(i) = specificenergy(norm(v(i,:)), mu_E, norm(r(i,:)), []);
end

hold on
plot((T./(t_sim).*365),specificEnergy);
hold off
xlabel('Time [days]'); ylabel('$\epsilon [km^2/s^2]$');
title('Specific Energy');
legend('Specific Energy');
grid on;


% Question 2

figure(2)
% Perturbed data
e_v = zeros(length(T), 3);
e_module = zeros(length(T), 1);
for i = 1:length(T)
    [e_v(i,:), e_module(i), ~] = eccentricity(r(i,:), v(i,:), h_v(i,:), mu_E);
end

hold on;
% Perturbed
plot(T,e_v(:,1), 'b-');
plot(T,e_v(:,2), 'r-');
plot(T,e_v(:,3), 'g-');
plot(T,e_module, 'k-', 'LineWidth', 2);

hold off;
xlabel('Time [days]'); ylabel('e [-]');
title('Eccentricity Vector');
legend('{$e_x$}, J2-perturbed','{$e_y$} J2-perturbed','{$z_p$}, J2-perturbed','${\| e \|$}, J2-perturbed', 'Interpreter', 'latex');
grid on;


%% Question 3
clc
clearvars
addpath("Functions 2bp\");

a = 7152; %km
e = 0.95;
num_cycles = 0.8; % Example input for the number of cycles (integer)

mu_E = astroConstants(13);
T_period = orbital_period(mu_E, a); % Orbital period [s]
t_sim = num_cycles * T_period; % Total simulation time
tspan = linspace(0, t_sim, 1000); % Time vector spanning num_cycles

x = zeros(length(tspan), length(e));

for j = 1:length(e)
    for i = 1:length(tspan)
        t = tspan(i);
        k = floor(t / T_period); % Number of full cycles
        t_remainder = rem(t, T_period); % Remainder time within the current cycle
        fun = @(E) keplertimelaw(mu_E, e(j), a, E, t_remainder);
        E = fsolve(fun, 0); % Solve for E using fsolve with initial guess 0
        x(i, j) = E + 2 * pi * k; % Total eccentric anomaly including full cycles
    end
end

% Convert eccentric anomaly to degrees for plotting
x_deg = rad2deg(x);

% Plot number of cycles vs eccentric anomaly
figure(3)

plot(tspan / T_period, x_deg);

xlabel('Number of Cycles');
ylabel('Eccentric Anomaly (deg)');
title('Eccentric Anomaly vs Number of Cycles');
legend('Eccentric Anomaly [°]');
grid on;