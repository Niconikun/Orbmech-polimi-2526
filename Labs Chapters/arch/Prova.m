addpath("Functions\");
%% Calculations

r = Y(:,1:3);
v = Y(:,4:6);

r_norm = zeros(length(t_vec),1);

% Normalization
for i = 1:length(t_vec)
    r_norm(i) = norm(r(i,:));
end

% Declination
delta = asin(r(:,3)./r_norm);

% Right ascension
for i = 1:length(t_vec)
    if (r(i,2)/r_norm(i) > 0)
        alpha = acos(r(i,1)/(r_norm*cos(delta(i))));
    elseif (r(i,2)/r_norm(i) <= 0)
        alpha = 2*pi - acos(r(i,1)/(r_norm*cos(delta(i))));
    end
end

% Latitude 
phi = delta;

% Longitude
theta_G = zeros(length(t_vec), 1);
for i = 1:length(t_vec)
    % angular velocity of earth = 7.2921*10^-5 rad/s
    theta_G(i) = (7.2921*10^-5)*t_vec(i);
end
lambda = alpha' - theta_G; 
lambda = rad2deg(lambda);

%% Plotting - Works with scatter plot (no background)

% Adjusting longitude for periodicity
lambda_periodic = mod(lambda + 180, 360) - 180; % Wraps longitude to [-180, 180]

% Creating a 2D matrix with adjusted longitude and corresponding latitude
ground_track = [lambda_periodic', phi'];

figure()
scatter(lambda_periodic, rad2deg(phi))
xlabel('Longitude [deg]'); ylabel('Latitude [deg]')
title('Spacecraft Ground Track');
grid on;
