clc
clear all
close all

%% Constants

deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;


addpath("ephemerides\")
addpath("Functions 2bp\")
addpath("time\")

AU = 149597870.691;              % Astronomical Unit [km]

% Sun
mu_Sun = astroConstants(4);      % Sun's gravitational parameter [km^3/s^2];
R_Sun = astroConstants(3);       % Sun's radius [km]

% Earth
mu_Earth = astroConstants(13);      % Sun's gravitational parameter [km^3/s^2];
R_Earth = astroConstants(23);       % Sun's radius [km]
v_Earth = 29.78; %km/s            % Earth orbital velocity around the Sun [km/s]
v_inf_Earth = 11.186; %km/s        % Earth escape velocity [km/s]

% Mars
mu_Mars = astroConstants(14);      % Sun's gravitational parameter [km^3/s^2];
R_Mars = astroConstants(24);       % Sun's radius [km]
v_Mars = 24.07; %km/s              % Mars orbital velocity around the Sun [km/s]
v_inf_Mars = 5.027; %km/s          % Mars escape velocity [km/s]


orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 0;  % Reduced to avoid warnings during computation


%% Time Window

%Earliest Departure
t_1_e_mjd2000 = date2mjd2000([2030,1,1,0,0,0]);

%---

t_1_l_mjd2000 = date2mjd2000([2003,8,1,0,0,0]);

% Reduced resolution for testing
n_points = 60;  % You can increase this later
departure_dates = linspace(t_1_e_mjd2000, t_1_l_mjd2000, n_points);

t_2_e_mjd2000 = date2mjd2000([2003,9,1,0,0,0]);
t_2_l_mjd2000 = date2mjd2000([2004,3,1,0,0,0]);

arrival_dates = linspace(t_2_e_mjd2000, t_2_l_mjd2000, n_points);


% Create meshgrid for contour plot
[Departure, Arrival] = meshgrid(departure_dates, arrival_dates);
delta_v = zeros(size(Departure));

fprintf('Computing porkchop plot...\n');

% Progress tracking
progress_interval = max(1, floor(n_points/10));


%% Lambert 1

for i = 1:n_points
    if mod(i, progress_interval) == 0
        fprintf('Progress: %.0f%%\n', (i/n_points)*100);
    end
    
    for j = 1:n_points
        % Get departure and arrival dates
        t_dep = departure_dates(i);
        t_arr = arrival_dates(j);
        
        % Calculate time of flight in days
        tof_days = t_arr - t_dep;
        
        % Skip if arrival is before departure or ToF is unrealistic
        %if tof_days <= 100 || tof_days > 600  % Reasonable bounds for Earth-Mars
        %    delta_v(j,i) = NaN;
        %    continue;
        %end
        
        tof_seconds = tof_days * 24 * 3600;
        
        try
            % Get planetary positions - FIXED: Ensure consistent coordinate systems
            [kep_earth, ~] = uplanet(t_dep, 3);  % Earth = 3
            [r1, v1] = kep2car(kep_earth, mu_Sun);
            
            [kep_mars, ~] = uplanet(t_arr, 4);   % Mars = 4
            [r2, v2] = kep2car(kep_mars, mu_Sun);
            
            % Ensure vectors are column vectors
            r1 = r1(:); v1 = v1(:);
            r2 = r2(:); v2 = v2(:);
            
            % Solve Lambert's problem
            [~, ~, ~, ERROR, VI, VF, ~, ~] = lambertMR(r1, r2, tof_seconds, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);
            
            if ERROR == 0
                % FIXED: Proper delta-v calculation with consistent vector dimensions
                cost_1 = norm(VI(:) - v1(:));
                cost_2 = norm(VF(:) - v2(:));
                delta_v(j,i) = cost_1 + cost_2;
            else
                delta_v(j,i) = NaN;
            end
            
        catch ME
            delta_v(j,i) = NaN;
        end
    end
end

%% FLyby (Hyperbolic Trajectory)


%% Lambert 2

