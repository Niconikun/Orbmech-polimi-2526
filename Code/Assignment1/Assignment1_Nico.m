%% Assignment1.m - Enhanced version with corrections and plotting
clc
clearvars
close all

%% Constants and Setup
deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;


% Add paths
addpath("ephemerides\")
addpath("Functions\")
addpath("Functions\timeConversion\")

SunImage = imread('img\sun.jpg');
EarthImage = imread('img\earth.jpg');
MarsImage = imread('img\mars.jpg');

AU = astroConstants(2);              % Astronomical Unit [km]

% Celestial bodies gravitational parameters and radii
mu_Sun = astroConstants(4);          % Sun's gravitational parameter [km^3/s^2]
R_Sun = astroConstants(3);           % Sun's radius [km]

mu_Earth = astroConstants(13);       % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23);        % Earth's radius [km]

mu_Mars = astroConstants(14);        % Mars's gravitational parameter [km^3/s^2]
R_Mars = astroConstants(24);         % Mars's radius [km]

% Lambert solver parameters
orbitType = 0;      % 0: direct transfer, 1: retrograde
Nrev = 0;           % Number of revolutions
Ncase = 0;          % 0: small-a option, 1: large-a option
optionsLMR = 0;     % Display options

%% Time Window Definition
% Reduced resolution for initial testing (increase later)
n_points = 120;

% Departure window from Earth (2030-2060)
t_dep_e = date2mjd2000([2046, 4, 20, 0, 0, 0]);
t_dep_l = date2mjd2000([2047, 5, 21, 0, 0, 0]);
departure_dates = linspace(t_dep_e, t_dep_l, n_points);

% Gravity assist at Mars
t_ga_e = date2mjd2000([2047, 4, 29, 0, 0, 0]);
t_ga_l = date2mjd2000([2048, 5, 30, 0, 0, 0]);
ga_dates = linspace(t_ga_e, t_ga_l, n_points);

% Arrival at asteroid
t_arr_e = date2mjd2000([2052, 6, 27, 0, 0, 0]);
t_arr_l = date2mjd2000([2053, 7, 27, 0, 0, 0]);
arrival_dates = linspace(t_arr_e, t_arr_l, n_points);

% Initialize delta-v matrix and other metrics
delta_v_total = zeros(n_points, n_points, n_points);
flyby_feasible = false(n_points, n_points, n_points);
min_periapsis = inf(n_points, n_points, n_points);

fprintf('Computing porkchop plot...\n');
progress_interval = max(1, floor(n_points^3/100));
counter = 0;

tof1days = zeros(n_points,n_points,n_points);
tof2days = zeros(n_points,n_points,n_points);

dv_flyby = zeros(n_points,n_points,n_points);
dv_departure = zeros(n_points,n_points,n_points);
dv_arrival = zeros(n_points,n_points,n_points);

%% Main Computation Loop
for i = 1:n_points      % Departure index
    for j = 1:n_points  % Gravity assist index
        for k = 1:n_points  % Arrival index
            
            counter = counter + 1;
            if mod(counter, progress_interval) == 0
                fprintf('Progress: %.1f%%\n', (counter/n_points^3)*100);
            end
            
            % Get dates
            t_dep = departure_dates(i);
            t_ga = ga_dates(j);
            t_arr = arrival_dates(k);

            tof1days(i,j,k) = t_ga - t_dep;
            tof2days(i,j,k) = t_arr - t_ga;
            
            % Check date order validity
            if ~(t_dep < t_ga && t_ga < t_arr)
                delta_v_total(i,j,k) = NaN;
                %fprintf("Invalid Date Combination\n");
                continue;
            end
            
            % Calculate time of flights
            tof_depga = (t_ga - t_dep) * 86400;  % seconds
            tof_gaarr = (t_arr - t_ga) * 86400;  % seconds
            
            try
                % Get planetary positions at respective times
                % Earth at departure
                kep_earth = uplanet(t_dep, 3);  % Earth ID = 3
                [r_earth, v_earth] = kep2car(kep_earth, mu_Sun);
                
                % Mars at gravity assist
                kep_mars = uplanet(t_ga, 4);    % Mars ID = 4
                [r_mars, v_mars] = kep2car(kep_mars, mu_Sun);
                
                % Asteroid at arrival (using example asteroid ID 257323)
                kep_aster = ephAsteroids(t_arr, 257323);
                [r_aster, v_aster] = kep2car(kep_aster, mu_Sun);
                
                % Ensure vectors are column vectors
                r_earth = r_earth(:);
                v_earth = v_earth(:);
                r_mars = r_mars(:);
                v_mars = v_mars(:);
                r_aster = r_aster(:);
                v_aster = v_aster(:);
                
                % Solve first Lambert arc (Earth to Mars)
                [~, ~, ~, error1, VI1, VF1, ~, ~] = lambertMR(r_earth, r_mars, tof_depga, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);
                
                % Solve second Lambert arc (Mars to Asteroid)
                [~, ~, ~, error2, VI2, VF2, ~, ~] = lambertMR(r_mars, r_aster, tof_gaarr, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);
                
                % Also ensure Lambert outputs are column vectors
                VI1 = VI1(:);
                VF1 = VF1(:);
                VI2 = VI2(:);
                VF2 = VF2(:);

                if error1 ~= 0 || error2 ~= 0
                    delta_v_total(i,j,k) = NaN;
                    %fprintf("No results found with this combination\n")
                    continue;
                end
                %fprintf("Combination found\n")
                % Flyby analysis at Mars
                v_inf_minus = VF1 - v_mars;    % Incoming hyperbolic excess velocity
                v_inf_plus = VI2 - v_mars;     % Outgoing hyperbolic excess velocity
                
                % Ensure vectors are column vectors with correct dimensions
                v_inf_minus = v_inf_minus(:);
                v_inf_plus = v_inf_plus(:);
                
                if length(v_inf_minus) ~= 3 || length(v_inf_plus) ~= 3
                    delta_v_total(i,j,k) = NaN;
                    %fprintf("Vector dimension error\n");
                    continue;
                end
                
                v_in_norm = norm(v_inf_minus);
                v_out_norm = norm(v_inf_plus);
                
                % Check for valid norms before calculating delta
                if v_in_norm == 0 || v_out_norm == 0
                    delta_v_total(i,j,k) = NaN;
                    %fprintf("Zero velocity norm\n");
                    continue;
                end
                
                % Calculate turning angle with safety check
                cos_delta = dot(v_inf_minus, v_inf_plus) / (v_in_norm * v_out_norm);
                % Ensure cos_delta is within valid range for acos
                cos_delta = max(-1, min(1, cos_delta));
                delta = acos(cos_delta);
                
                % Check if delta is a valid scalar
                if ~isscalar(delta) || ~isreal(delta)
                    delta_v_total(i,j,k) = NaN;
                    %fprintf("Invalid delta calculation\n");
                    continue;
                end
                
                % Solve for periapsis radius using fzero
                if delta > 0 && delta < pi
                    fun = @(r_p) asin(1/(1 + r_p*v_out_norm^2/mu_Mars)) + ...
                                 asin(1/(1 + r_p*v_in_norm^2/mu_Mars)) - delta;
                    
                    % Better initial guess based on turning angle
                    if delta < 30*pi/180
                        x0 = 5 * R_Mars;  % Small turning angle -> higher flyby
                    else
                        x0 = 2 * R_Mars;  % Larger turning angle -> lower flyby
                    end
                    
                    % Flyby safety parameters
                    safety_margin = 200;  % km - minimum altitude above Mars surface
                    max_flyby_altitude = 10000;  % km - maximum reasonable flyby altitude

                    try
                        options = optimset('TolX', 1e-12, 'TolFun', 1e-16, ...
                                          'MaxIter', 100, 'Display', 'off');
                        r_p = fzero(fun, x0, options);
                        
                        % Store minimum periapsis
                        min_periapsis(i,j,k) = r_p;
                        
                     
                        % Check flyby feasibility
                        % Check flyby feasibility with proper validation
                        if r_p > 0 && isfinite(r_p) && r_p >= R_Mars + safety_margin
                            flyby_feasible(i,j,k) = true;
                            
                            % Calculate velocities at periapsis
                            v_p_in = sqrt(v_in_norm^2 + 2*mu_Mars/r_p);
                            v_p_out = sqrt(v_out_norm^2 + 2*mu_Mars/r_p);
                            dv_flyby(i,j,k) = abs(v_p_out - v_p_in);
                            
                            % Calculate total delta-v
                            dv_departure(i,j,k) = norm(VI1 - v_earth);
                            dv_arrival(i,j,k) = norm(VF2 - v_aster);
                            delta_v_total(i,j,k) = dv_departure + dv_flyby + dv_arrival;

                            
                        else
                            flyby_feasible(i,j,k) = false;
                            delta_v_total(i,j,k) = NaN;
                        end
                    catch
                        delta_v_total(i,j,k) = NaN;
                    end
                else
                    delta_v_total(i,j,k) = NaN;
                end
                
            catch ME
                delta_v_total(i,j,k) = NaN;
                fprintf('Error at (%d,%d,%d): %s\n', i, j, k, ME.message);
            end
        end
    end
end

fprintf('Computation complete!\n');

%% Find Optimal Trajectory
% Find minimum delta-v among feasible trajectories
feasible_indices = find(flyby_feasible);
if ~isempty(feasible_indices)
    [min_dv, min_idx] = min(delta_v_total(feasible_indices));
    [i_opt, j_opt, k_opt] = ind2sub(size(delta_v_total), feasible_indices(min_idx));
    
    % Optimal dates
    t_dep_opt = departure_dates(i_opt);
    t_ga_opt = ga_dates(j_opt);
    t_arr_opt = arrival_dates(k_opt);
    
    % Convert MJD2000 to date vectors using your function
    dep_date_vec = mjd20002date(t_dep_opt);  % Returns [year, month, day, hour, minute, second]
    ga_date_vec = mjd20002date(t_ga_opt);
    arr_date_vec = mjd20002date(t_arr_opt);
    
    % Create datetime objects
    dep_datetime = datetime(dep_date_vec(1), dep_date_vec(2), dep_date_vec(3), ...
                           dep_date_vec(4), dep_date_vec(5), floor(dep_date_vec(6)));
    ga_datetime = datetime(ga_date_vec(1), ga_date_vec(2), ga_date_vec(3), ...
                          ga_date_vec(4), ga_date_vec(5), floor(ga_date_vec(6)));
    arr_datetime = datetime(arr_date_vec(1), arr_date_vec(2), arr_date_vec(3), ...
                           arr_date_vec(4), arr_date_vec(5), floor(arr_date_vec(6)));
    
    fprintf('\n=== OPTIMAL TRAJECTORY FOUND ===\n');
    fprintf('Departure: %s (MJD2000 = %.2f)\n', datestr(dep_datetime), t_dep_opt);
    fprintf('Gravity Assist: %s (MJD2000 = %.2f)\n', datestr(ga_datetime), t_ga_opt);
    fprintf('Arrival: %s (MJD2000 = %.2f)\n', datestr(arr_datetime), t_arr_opt);
    fprintf('Total ΔV: %.4f km/s\n', min_dv);
    fprintf('ΔV Departure: %.4f km/s\n', dv_departure(i_opt,j_opt,k_opt));
    fprintf('ΔV Powered Gravity Assist: %.4f km/s\n', dv_flyby(i_opt,j_opt,k_opt));
    fprintf('ΔV Arrival: %.4f km/s\n', dv_arrival(i_opt,j_opt,k_opt));
    fprintf('Time of Flight: %.2f days\n', t_arr_opt - t_dep_opt);
    fprintf('Flyby periapsis: %.2f km (%.2f R_Mars)\n', ...
            min_periapsis(i_opt,j_opt,k_opt), ...
            min_periapsis(i_opt,j_opt,k_opt)/R_Mars);
    
    % Display additional trajectory information
    fprintf('\n=== TRAJECTORY DETAILS ===\n');
    fprintf('Earth-Mars transfer time: %.2f days\n', t_ga_opt - t_dep_opt);
    fprintf('Mars-Asteroid transfer time: %.2f days\n', t_arr_opt - t_ga_opt);
    fprintf('Flyby altitude: %.2f km\n', min_periapsis(i_opt,j_opt,k_opt) - R_Mars);        
    
    % Times of flight
    tof_depga = (t_ga_opt - t_dep_opt) * 86400;
    tof_gaarr = (t_arr_opt - t_ga_opt) * 86400;
    
    % Get positions and velocities
    kep_earth = uplanet(t_dep_opt, 3);
    [r_earth, v_earth] = kep2car(kep_earth, mu_Sun);
    
    kep_mars = uplanet(t_ga_opt, 4);
    [r_mars, v_mars] = kep2car(kep_mars, mu_Sun);
    
    kep_aster = ephAsteroids(t_arr_opt, 257323);
    [r_aster, v_aster] = kep2car(kep_aster, mu_Sun);

    % Solve Lambert problems
    [~, ~, ~, ~, VI1, VF1, ~, ~] = lambertMR(r_earth, r_mars, tof_depga, mu_Sun, 0, 0, 0, 0);
    [~, ~, ~, ~, VI2, VF2, ~, ~] = lambertMR(r_mars, r_aster, tof_gaarr, mu_Sun, 0, 0, 0, 0);
    
    % Flyby calculations
    v_inf_minus = VF1 - v_mars';
    v_inf_plus = VI2 - v_mars';
    
    % Ensure column vectors
    v_inf_minus = v_inf_minus(:);
    v_inf_plus = v_inf_plus(:);
    
    v_in_norm = norm(v_inf_minus);
    v_out_norm = norm(v_inf_plus);
    
    delta = acos(dot(v_inf_minus, v_inf_plus) / (v_in_norm * v_out_norm));

    opt_trajectory.t_dep = t_dep_opt;
    opt_trajectory.t_ga = t_ga_opt;
    opt_trajectory.t_arr = t_arr_opt;
    opt_trajectory.r_earth = r_earth;
    opt_trajectory.r_mars = r_mars;
    opt_trajectory.r_aster = r_aster;
    opt_trajectory.VI1 = VI1;
    opt_trajectory.VF1 = VF1;
    opt_trajectory.VI2 = VI2;
    opt_trajectory.VF2 = VF2;
    opt_trajectory.v_earth = v_earth;
    opt_trajectory.v_mars = v_mars;
    opt_trajectory.v_aster = v_aster;
    
    % Store flyby data
    opt_flyby.r_p = min_periapsis(i_opt,j_opt,k_opt);
    opt_flyby.v_inf_minus = v_inf_minus;
    opt_flyby.v_inf_plus = v_inf_plus;
    opt_flyby.v_in_norm = v_in_norm;
    opt_flyby.v_out_norm = v_out_norm;
    opt_flyby.delta = delta;
    opt_flyby.r_mars = r_mars;
    opt_flyby.v_mars = v_mars;
    
    %[opt_trajectory, opt_flyby] = compute_optimal_trajectory(t_dep_opt, t_ga_opt, t_arr_opt);
    
else
    fprintf('\nNo feasible trajectories found!\n');
    return;
end

%% Plotting Section
fprintf('\nGenerating plots...\n');

% 1. Porkchop Plot - Departure vs Gravity Assist (for fixed arrival)
figure('Name', 'Porkchop Plot');

% Subplot 1: Departure vs Gravity Assist
subplot(1, 3, 1);
% Extract 2D slice for fixed arrival index

dv_slice_dep_ga = squeeze(dv_departure(:,:,k_opt));

create_porkchop_plot(departure_dates, ga_dates, dv_slice_dep_ga, i_opt, j_opt, tof1days(:,:,k_opt), dv_departure(i_opt,j_opt,k_opt));
title('Departure vs Gravity Assist (Fixed Arrival)');

% Subplot 2: Gravity Assist vs Arrival
subplot(1, 3, 2);
% Extract 2D slice for fixed departure index

dv_slice_ga_arr = squeeze(dv_arrival(i_opt,:,:));

create_porkchop_plot(ga_dates, arrival_dates, dv_slice_ga_arr, j_opt, k_opt, tof2days(i_opt,:,:), dv_arrival(i_opt,j_opt,k_opt));
title('Gravity Assist vs Arrival (Fixed Departure)');

% Subplot 3: Departure vs Arrival
subplot(1, 3, 3);
% Extract 2D slice for fixed departure index

dv_slice_dep_arr = squeeze(delta_v_total(:,j_opt,:));

create_porkchop_plot(departure_dates, arrival_dates, dv_slice_dep_arr, i_opt, k_opt, tof2days(:,j_opt,:), delta_v_total(i_opt,j_opt,k_opt));
title('Gravity Assist vs Arrival (Fixed Departure)');

%% New Figure: Comprehensive Trajectory Visualization
fprintf('\nGenerating comprehensive trajectory plots...\n');

figure('Name', 'Optimized Trajectory Analysis');

% ------------------------------------------
% Subplot 1: 3D view of the optimized trajectory
subplot(2, 3, [1,5]);
plot_3d_trajectory_comprehensive(opt_trajectory, opt_flyby,n_points,SunImage);
title('3D View: Complete Transfer Orbit');

% ------------------------------------------
% Subplot 2: 2D projection (ecliptic plane)
subplot(1, 3, 3);
plot_2d_projection(opt_trajectory, opt_flyby,n_points);
title('2D Projection: Ecliptic Plane View');

% Add an overall title for the figure
sgtitle('Optimized Trajectory: Earth → Mars → Asteroid', 'FontSize', 14, 'FontWeight', 'bold');

% ------------------------------------------
% Plot 3: 3D view of Mars flyby

figure('Name','Flyby');
plot_3d_flyby(opt_flyby,MarsImage,v_mars);
title('3D View: Mars Flyby Geometry');


%% Helper Functions for the New Plots

function plot_3d_trajectory_comprehensive(traj, ~,~,img1)
    % Comprehensive 3D view of the complete trajectory
    
    % Constants

    mu_Sun = astroConstants(4);
    R_Sun = astroConstants(3);
    
    hold on;
    grid on;
    box on;
    
    scaleFactor = 40;

    % Set 3D viewing properties
    view(45, 30);  % 3D perspective
    axis equal;
    
    % Plot the Sun at origin
    [X_sun, Y_sun, Z_sun] = sphere(50);
    surf(R_Sun * X_sun*scaleFactor, R_Sun * Y_sun*scaleFactor, R_Sun * Z_sun*scaleFactor, 'FaceColor', 'texturemap', 'CData', img1, 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName','Sun');
    
    % Earth orbit
    
    [kep_earth] = car2kep(traj.r_earth, traj.v_earth, mu_Sun);
    
    a_earth = kep_earth(1);

    y_earth = [traj.r_earth; traj.v_earth];
    T_earth = 2*pi*sqrt(a_earth^3/mu_Sun); % Orbital period [s]
    tspan_earth = linspace( 0, T_earth, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Yearth ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_earth, y_earth, options );
    
    % Calculate ground track
    earth_orbit = Yearth(:,1:3);

    plot3(earth_orbit(:,1), earth_orbit(:,2), earth_orbit(:,3), ...
          'b--', 'LineWidth', 0.5, 'Color', [0.3, 0.3, 1],'DisplayName','Earth Orbit');
    
    % Mars orbit (circular approximation)
    [kep_mars] = car2kep(traj.r_mars, traj.v_mars, mu_Sun);
    a_mars = kep_mars(1);
    
    y_mars = [traj.r_mars; traj.v_mars];
    T_mars = 2*pi*sqrt(a_mars^3/mu_Sun); % Orbital period [s]
    tspan_mars = linspace( 0, T_mars, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ymars ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_mars, y_mars, options );
    
    % Calculate ground track
    mars_orbit = Ymars(:,1:3);

    plot3(mars_orbit(:,1), mars_orbit(:,2), mars_orbit(:,3), ...
          'r--', 'LineWidth', 0.5, 'Color', [1, 0.3, 0.3],'DisplayName','Mars Orbit');

     % Asteroid orbit
    
    [kep_aster] = car2kep(traj.r_aster, traj.v_aster, mu_Sun);
    
    a_aster = kep_aster(1);

    y_aster = [traj.r_aster; traj.v_aster];
    T_aster = 2*pi*sqrt(a_aster^3/mu_Sun); % Orbital period [s]
    tspan_aster = linspace( 0, T_aster, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Yaster ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_aster, y_aster, options );
    
    % Calculate ground track
    aster_orbit = Yaster(:,1:3);

    plot3(aster_orbit(:,1), aster_orbit(:,2), aster_orbit(:,3), ...
          'b--', 'LineWidth', 0.5, 'Color', [0.3, 1, 0.3], 'DisplayName','Asteroid Orbit');

    % Plot planetary positions
    % Earth
    plot3(traj.r_earth(1), traj.r_earth(2), traj.r_earth(3), ...
          'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
          'LineWidth', 2, 'DisplayName','Earth');
    
    % Mars
    plot3(traj.r_mars(1), traj.r_mars(2), traj.r_mars(3), ...
          'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
          'LineWidth', 2, 'DisplayName','Mars');
    
    % Asteroid
    plot3(traj.r_aster(1), traj.r_aster(2), traj.r_aster(3), ...
          'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', ...
          'LineWidth', 2,'DisplayName','Asteroid');
    
    % Earth to Mars arc (using ODE solver)
    y_transfer1 = [traj.r_earth; traj.VI1'];
    T_transfer1 = (traj.t_ga - traj.t_dep)*24*60*60; % Orbital period [s]
    
    tspan_transfer1 = linspace(0, T_transfer1, 1000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ytransfer1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_transfer1, y_transfer1, options );
            
    r_arc1 = Ytransfer1(:,1:3);

    
    % Mars to Asteroid arc
    y_transfer2 = [traj.r_mars; traj.VI2'];
    T_transfer2 = (traj.t_arr - traj.t_ga)*24*60*60; % Orbital period [s]
    
    tspan_transfer2 = linspace(0, T_transfer2, 1000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ytransfer2 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_transfer2, y_transfer2, options );
            
    r_arc2 = Ytransfer2(:,1:3);

    plot3(r_arc1(:,1), r_arc1(:,2), r_arc1(:,3), 'b-', 'LineWidth', 1,'DisplayName','Earth-Mars Arc');
    plot3(r_arc2(:,1), r_arc2(:,2), r_arc2(:,3), 'r-', 'LineWidth', 1,'DisplayName','Mars-Asteroid Arc');
    
    % Add labels
    text(traj.r_earth(1), traj.r_earth(2), traj.r_earth(3), ...
         ' Earth', 'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
    text(traj.r_mars(1), traj.r_mars(2), traj.r_mars(3), ...
         ' Mars', 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    text(traj.r_aster(1), traj.r_aster(2), traj.r_aster(3), ...
         ' Asteroid', 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
    text(0, 0, 0, ' Sun', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
    
    % Add legend (simplified)
    legend('Location', 'best', 'FontSize', 7);
    
    % Labels and formatting
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    
    % Set appropriate axis limits
    max_dist = max([norm(traj.r_earth), norm(traj.r_mars), norm(traj.r_aster)]) * 1.1;
    xlim([-max_dist, max_dist]);
    ylim([-max_dist, max_dist]);
    zlim([-max_dist/5, max_dist/5]);  % Flatten for better view
    
    % Add lighting
    camlight('headlight');
    lighting gouraud;
    
    hold off;
end

function plot_2d_projection(traj, ~,~)
    % 2D projection (ecliptic plane view)
    
    mu_Sun = astroConstants(4);
    
    hold on;
    grid on;
    box on;
    
    % Plot the Sun at origin
    plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'y','DisplayName','Sun');
    
    % Generate transfer arcs

    
    % Earth to Mars arc (using ODE solver)
    y_transfer1 = [traj.r_earth; traj.VI1'];
    T_transfer1 = (traj.t_ga - traj.t_dep)*24*60*60; % Orbital period [s]
    
    tspan_transfer1 = linspace(0, T_transfer1, 1000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ytransfer1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_transfer1, y_transfer1, options );
            
    r_arc1 = Ytransfer1(:,1:3);

    
    % Mars to Asteroid arc
    y_transfer2 = [traj.r_mars; traj.VI2'];
    T_transfer2 = (traj.t_arr - traj.t_ga)*24*60*60; % Orbital period [s]
    
    tspan_transfer2 = linspace(0, T_transfer2, 1000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ytransfer2 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_transfer2, y_transfer2, options );
            
    r_arc2 = Ytransfer2(:,1:3);
    
     % Earth orbit
    
    [kep_earth] = car2kep(traj.r_earth, traj.v_earth, mu_Sun);
    
    a_earth = kep_earth(1);

    y_earth = [traj.r_earth; traj.v_earth];
    T_earth = 2*pi*sqrt(a_earth^3/mu_Sun); % Orbital period [s]
    tspan_earth = linspace( 0, T_earth, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Yearth ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_earth, y_earth, options );
    
    % Calculate ground track
    earth_orbit = Yearth(:,1:3);

    plot(earth_orbit(:,1), earth_orbit(:,2), ...
          'b--', 'LineWidth', 0.5, 'Color', [0.3, 0.3, 1],'DisplayName','Earth Orbit');
    
    % Mars orbit
    [kep_mars] = car2kep(traj.r_mars, traj.v_mars, mu_Sun);
    a_mars = kep_mars(1);
    
    y_mars = [traj.r_mars; traj.v_mars];
    T_mars = 2*pi*sqrt(a_mars^3/mu_Sun); % Orbital period [s]
    tspan_mars = linspace( 0, T_mars, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Ymars ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_mars, y_mars, options );
    
    % Calculate ground track
    mars_orbit = Ymars(:,1:3);

    plot(mars_orbit(:,1), mars_orbit(:,2), ...
          'r--', 'LineWidth', 0.5, 'Color', [1, 0.3, 0.3],'DisplayName','Mars Orbit');

     % Asteroid orbit
    
    [kep_aster] = car2kep(traj.r_aster, traj.v_aster, mu_Sun);
    
    a_aster = kep_aster(1);

    y_aster = [traj.r_aster; traj.v_aster];
    T_aster = 2*pi*sqrt(a_aster^3/mu_Sun); % Orbital period [s]
    tspan_aster = linspace( 0, T_aster, 10000); % Reduced to 10 periods for clarity
    
    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-11 );
    
    % Perform the integration
    [ ~, Yaster ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_aster, y_aster, options );
    
    % Calculate ground track
    aster_orbit = Yaster(:,1:3);

    plot(aster_orbit(:,1), aster_orbit(:,2), ...
          'b--', 'LineWidth', 0.5, 'Color', [0.3, 1, 0.3],'DisplayName','Asteroid Orbit');
    
    % Plot transfer arcs
    plot(r_arc1(:,1), r_arc1(:,2), 'b-', 'LineWidth', 2, 'DisplayName','Earth-Mars Arc');
    plot(r_arc2(:,1), r_arc2(:,2), 'r-', 'LineWidth', 2, 'DisplayName','Mars-Asteroid Arc');
    
    % Plot planetary positions
    plot(traj.r_earth(1), traj.r_earth(2), 'bo', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'b', 'LineWidth', 2, 'DisplayName','Earth');
    plot(traj.r_mars(1), traj.r_mars(2), 'ro', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName','Mars');
    plot(traj.r_aster(1), traj.r_aster(2), 'go', ...
         'MarkerSize', 6, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'DisplayName','Asteroid');
    
    
    % Add labels
    text(traj.r_earth(1), traj.r_earth(2), ' Earth', ...
         'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
    text(traj.r_mars(1), traj.r_mars(2), ' Mars', ...
         'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    text(traj.r_aster(1), traj.r_aster(2), ' Asteroid', ...
         'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
    text(0, 0, ' Sun', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
    
    % Formatting
    axis equal;
    xlabel('X (km)');
    ylabel('Y (km)');
    
    % Set axis limits
    max_dist = max([a_earth, a_mars, norm(traj.r_aster)]) * 1.1;
    xlim([-max_dist, max_dist]);
    ylim([-max_dist, max_dist]);
    
    % Add grid
    grid on;
    grid minor;
    
    % Add legend
    legend('Location', 'best', 'FontSize', 7);
    
    hold off;
end

function plot_3d_flyby(flyby, img, V_T)
    % 3D view of the Mars flyby geometry
    
    mu_Mars = astroConstants(14);
    R_Mars = astroConstants(24);
    
    % Extract flyby parameters
    v_inf_minus = flyby.v_inf_minus;
    v_inf_plus = flyby.v_inf_plus;
    v_in_norm = flyby.v_in_norm;
    v_out_norm = flyby.v_out_norm;
    r_p = flyby.r_p;
    
    % Calculate hyperbola parameters
    a_in  = -mu_Mars / v_in_norm^2;
    a_out = -mu_Mars / v_out_norm^2;
    
    e_in  = 1 - r_p / a_in;
    e_out = 1 - r_p / a_out;
    
    delta_in  = 2 * asin(1 / e_in);
    delta_out = 2 * asin(1 / e_out);
    
    % Get orbital elements using fixed plane function
    h_dir = cross(v_inf_minus, v_inf_plus) / norm(cross(v_inf_minus, v_inf_plus));
    [i_in, Omega_in, omega_in] = plane_elements(v_inf_minus, v_inf_plus, delta_in, h_dir);
    [i_out, Omega_out, omega_out] = plane_elements(v_inf_minus, v_inf_plus, delta_out, h_dir);
    
    % Create Keplerian elements
    kep_in  = [a_in, e_in, i_in, Omega_in, omega_in, 0];
    kep_out = [a_out, e_out, i_out, Omega_out, omega_out, 0];

    c_in  = abs( a_in  * e_in );
    c_out = abs( a_out * e_out );

    h_dir = cross(v_inf_minus,v_inf_plus) / norm( cross(v_inf_minus,v_inf_plus) );

    r_p_dir = -rodrigues_rotation(h_dir,v_inf_plus,pi/2-delta_out/2) / v_out_norm;

    r_p_in  = c_in  * r_p_dir;
    r_p_out = c_out * r_p_dir;

    asymptote_in  = r_p_in  + v_inf_minus * [-1 -1000];
    asymptote_out = r_p_out + v_inf_plus  * [ 1  1000];

    tratteggio = r_p_in * [1 -2.5];
    
    hold on
    
    % Plot Mars sphere
    [X_mars, Y_mars, Z_mars] = sphere(50);
    surf(R_Mars * X_mars, R_Mars * Y_mars, R_Mars * Z_mars, ...
         'FaceColor', 'texturemap', 'CData', img, 'EdgeColor', 'none');
    
    my_plotOrbit(kep_in , mu_Mars, 500, -1.4,    0, 'b')
    my_plotOrbit(kep_out, mu_Mars, 500, 0   , +1.4, 'r')
    plot3(asymptote_in(1,:),asymptote_in(2,:),asymptote_in(3,:),'LineWidth',1.5)
    plot3(asymptote_out(1,:),asymptote_out(2,:),asymptote_out(3,:),'LineWidth',1.5)
    plot3(tratteggio(1,:),tratteggio(2,:),tratteggio(3,:),'--k','LineWidth',1.5)
    my_plotOrbit(kep_out, mu_Mars, 1   , 0   ,    0, 'ob')
    quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
    axis equal
    grid on
    box on;
    legend('Incoming hyperbola', 'Outcoming hyperbola','Asymptote of incoming hyperbola','Asymptote of outcoming hyperbola')
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    view(3);
    
    % Add lighting
    camlight('headlight');
    lighting gouraud;
    
    hold off;
end




%% Helper Functions for Plotting

function create_porkchop_plot(x_dates, y_dates, dv_slice, opt_x, opt_y,TOFdays, dv_min)
    % Create 2D porkchop plot
    
    TOFdays = squeeze(TOFdays);
    
    % Create meshgrid
    [X, Y] = meshgrid(x_dates, y_dates);
    
    % Create contour plot
    contourf(X, Y, dv_slice, 20, 'LineStyle', 'none');
    colormap jet;
    c = colorbar;
    %clim([min(dv_slice(:)) max(dv_slice(:))]);                 % *<- key line*
    c.Ticks = 5:0.5:10;
    c.Label.String = '\Delta v [km/s]';
    
    hold on; grid on;
    colorbar;
    colormap jet;
    xlabel('Departure date');
    ylabel('Arrival date');
    
    % Add contour lines
    ToF_levels = 50:500:5000;
    [C2,h2] = contour(X, Y, TOFdays, ToF_levels, 'k--', 'LineWidth', 0.8);
    clabel(C2,h2,'Color','k','FontSize',8);
    
    % Mark optimal point if indices are valid
    plot(x_dates(opt_x), y_dates(opt_y), 'wo', 'MarkerFaceColor','k', 'MarkerSize',7);
    text(x_dates(opt_x), y_dates(opt_y), sprintf('  min \\DeltaV = %.2f km/s', dv_min), ...
    'Color','w','FontSize',8,'FontWeight','bold');
    
    
    % Format dates on axes
    datetick('x','dd-mmm-yyyy','keeplimits');
    datetick('y','dd-mmm-yyyy','keeplimits');
    
    hold off;
end


%% Additional Helper Functions

function v_rot = rodrigues_rotation(u,v,delta)
    % Rodrigues rotation formula
    v_rot = v*cos(delta) + cross(u,v)*sin(delta) + u*dot(u,v)*(1-cos(delta));
end


function [i, Omega, omega] = plane_elements(~, v_inf_plus, delta, h_dir)
    % Calculate orbital plane parameters from flyby
    
    tol = 1e-12;
    x_vec = [1; 0; 0];
    z_vec = [0; 0; 1];
    
    % Calculate inclination
    i = acos(h_dir(3));
    
    % Calculate longitude of ascending node (Omega)
    N = cross(z_vec, h_dir);
    N_norm = norm(N);
    
    if N_norm < tol
        N = x_vec;
        N_norm = 1;
        Omega = 0;
    else
        if N(2) >= 0
            Omega = acos(N(1) / N_norm);
        else
            Omega = 2*pi - acos(N(1) / N_norm);
        end
    end
    
    % Calculate argument of periapsis (omega)
    beta = (pi - delta) / 2;
    v_inf_plus_norm = v_inf_plus / norm(v_inf_plus);
    e_vec_norm = -rodrigues_rotation(h_dir, v_inf_plus_norm, beta);
    
    if (i > tol) && (i < pi - tol)
        if e_vec_norm(3) >= 0
            omega = acos(dot(N, e_vec_norm) / N_norm);
        else
            omega = 2*pi - acos(dot(N, e_vec_norm) / N_norm);
        end
    elseif i <= tol
        if e_vec_norm(2) >= 0
            omega = acos(dot(N, e_vec_norm) / N_norm);
        else
            omega = 2*pi - acos(dot(N, e_vec_norm) / N_norm);
        end
    else  % i >= pi - tol
        if e_vec_norm(2) <= 0
            omega = acos(dot(N, e_vec_norm) / N_norm);
        else
            omega = 2*pi - acos(dot(N, e_vec_norm) / N_norm);
        end
    end
end



function my_plotOrbit(kep, mu, npoints, theta_i, theta_f, string)

theta = linspace(theta_i,theta_f,npoints);

R = zeros(3,npoints);

for i = 1:npoints

    kep(6) = theta(i);

    r = kep2car(kep,mu);

    R(:,i) = r;

end

plot3(R(1,:),R(2,:),R(3,:),string,'LineWidth',1.5)

end
