%% fsolve

clc
clearvars
addpath("Functions 2bp\");

a = 7000;
e = [0 0.2 0.4 0.6 0.8 0.95];

mu_E = astroConstants(13);
num_cycles = 1; % Example input for the number of cycles (integer)
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
figure(1)

plot(tspan / T_period, x_deg, 'DisplayName', ['e = ' num2str(e(j))]);

xlabel('Number of Cycles');
ylabel('Eccentric Anomaly (deg)');
title('Eccentric Anomaly vs Number of Cycles');
legend('Location', 'best');
grid on;

figure(2)
[T, E] = meshgrid(tspan / T_period, e);
surf(T, E, x_deg');
xlabel('Number of Cycles');
ylabel('Eccentricity');
zlabel('Eccentric Anomaly (deg)');
title('Eccentric Anomaly vs Number of Cycles and Eccentricity');
colorbar;
