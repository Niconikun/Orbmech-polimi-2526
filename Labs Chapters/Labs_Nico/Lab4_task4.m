clc
clearvars
close all

deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;

N_planet = 1; % Number of planet for uplanet

switch N_planet %in km/s
    case 1
        delta_v1_max = 7;
        planetName = "Mercury";
        departure_time_earliest = [2023,11,1,0,0,0];
        departure_time_latest = [2025,1,1,0,0,0];
        arrival_time_earliest = [2024,4,1,0,0,0];
        arrival_time_latest = [2025,3,1,0,0,0];
        R_planet = astroConstants(21);
        
    case 2
        delta_v1_max = 3;
        planetName = "Venus";
        departure_time_earliest = [2024,6,1,0,0,0];
        departure_time_latest = [2026,11,1,0,0,0];
        arrival_time_earliest = [2024,12,1,0,0,0];
        arrival_time_latest = [2027,6,1,0,0,0];
        R_planet = astroConstants(22);
    case 4
        delta_v1_max = 3.5;
        planetName = "Mars";
        departure_time_earliest = [2025,8,1,0,0,0];
        departure_time_latest = [2031,1,1,0,0,0];
        arrival_time_earliest =[2026,1,1,0,0,0];
        arrival_time_latest =[2032,1,1,0,0,0];
        R_planet = astroConstants(24);
    case 5
        delta_v1_max = 9.1;
        planetName = "Jupiter";
        departure_time_earliest =[2026,1,1,0,0,0];
        departure_time_latest = [2028,1,1,0,0,0];
        arrival_time_earliest =[2028,6,1,0,0,0];
        arrival_time_latest =[2034,1,1,0,0,0];
        R_planet = astroConstants(25);
    case 6
        delta_v1_max = 11.5;
        planetName = "Saturn";
        departure_time_earliest =[2027,9,1,0,0,0];
        departure_time_latest = [2029,10,1,0,0,0];
        arrival_time_earliest =[2030,4,1,0,0,0];
        arrival_time_latest =[2036,3,1,0,0,0];
        R_planet = astroConstants(26);
    case 7
        delta_v1_max = 12.1;
        planetName = "Uranus";
        departure_time_earliest =[2027,1,1,0,0,0];
        departure_time_latest = [2029,1,1,0,0,0];
        arrival_time_earliest =[2031,4,1,0,0,0];
        arrival_time_latest =[2045,12,1,0,0,0];
        R_planet = astroConstants(27);
    case 8
        delta_v1_max = 12.5;
        planetName = "Neptune";
        departure_time_earliest =[2025,1,1,0,0,0];
        departure_time_latest = [2026,10,1,0,0,0];
        arrival_time_earliest =[2036,1,1,0,0,0];
        arrival_time_latest =[2055,06,1,0,0,0];
        R_planet = astroConstants(28);
    case 9
        delta_v1_max = 12.5; %leave the same as neptune for fun
        planetName = "Pluto";
        departure_time_earliest =[2025,6,1,0,0,0];
        departure_time_latest = [2026,6,1,0,0,0];
        arrival_time_earliest =[2046,1,1,0,0,0];
        arrival_time_latest =[2125,6,1,0,0,0];
        R_planet = astroConstants(29);
end


addpath("ephemerides\")
addpath("Functions 2bp\")
addpath("time\")

mu_Sun = astroConstants(4);      % Sun's gravitational parameter [km^3/s^2];
R_Sun = astroConstants(3);       % Sun's radius [km]
AU = 149597870.691;              % Astronomical Unit [km]

R_E = astroConstants(23);
SunImage = imread('SunTexture.jpg');
EarthImage = imread('EarthTexture.jpg');
textureImage = imread(planetName+'Texture.jpg');
orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 0;  % Reduced to avoid warnings during computation

%% Time Window
t_1_e_mjd2000 = date2mjd2000(departure_time_earliest);
t_1_l_mjd2000 = date2mjd2000(departure_time_latest);

% Reduced resolution for testing
n_points = 800;  % You can increase this later
departure_dates = linspace(t_1_e_mjd2000, t_1_l_mjd2000, n_points);

t_2_e_mjd2000 = date2mjd2000(arrival_time_earliest);
t_2_l_mjd2000 = date2mjd2000(arrival_time_latest);

arrival_dates = linspace(t_2_e_mjd2000, t_2_l_mjd2000, n_points);

% Create meshgrid for contour plot
[Departure, Arrival] = meshgrid(departure_dates, arrival_dates);
delta_v = zeros(size(Departure));

fprintf('Computing porkchop plot...\n');

% Progress tracking
progress_interval = max(1, floor(n_points/10));

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
  
        
        tof_seconds = tof_days * 24 * 3600;
        
        try
            % Get planetary positions - FIXED: Ensure consistent coordinate systems
            [kep_earth, ~] = uplanet(t_dep, 3);
            kep_earth_deg = [kep_earth(1),kep_earth(2),rad2deg(kep_earth(3)),rad2deg(kep_earth(4)),rad2deg(kep_earth(5)),rad2deg(kep_earth(6))];
            [r1, v1] = kep2car(kep_earth_deg, mu_Sun);
            
            [kep_mars, ~] = uplanet(t_arr, N_planet);
            kep_mars_deg = [kep_mars(1),kep_mars(2),rad2deg(kep_mars(3)),rad2deg(kep_mars(4)),rad2deg(kep_mars(5)),rad2deg(kep_mars(6))];
            [r2, v2] = kep2car(kep_mars_deg, mu_Sun);
            
            % Ensure vectors are column vectors
            r1 = r1(:); v1 = v1(:);
            r2 = r2(:); v2 = v2(:);
            
            % Solve Lambert's problem
            [~, ~, ~, ERROR, VI, VF, ~, ~] = lambertMR(r1, r2, tof_seconds, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);

            VI = VI(:);
            VF = VF(:);
            
            if ERROR == 0
                % FIXED: Proper delta-v calculation with consistent vector dimensions
                cost_1 = norm(VI(:) - v1(:));
                %if cost_1 <= delta_v1_max
                    delta_v_departure = cost_1;
                %else
                %    delta_v_departure = delta_v1_max; % Set to max if exceeded
                %end
                delta_v_arrival = norm(VF(:) - v2(:));
                delta_v(j,i) = delta_v_departure + delta_v_arrival;
            else
                delta_v(j,i) = NaN;
            end
            
        catch ME
            delta_v(j,i) = NaN;
        end
    end
end

%% Create Porkchop Plot
fprintf('Creating porkchop plot...\n');

figure(1)

% Convert mjd2000 to datenum for plotting
departure_datenum = departure_dates + 730486.5;
arrival_datenum = arrival_dates + 730486.5;
[Departure_dn, Arrival_dn] = meshgrid(departure_datenum, arrival_datenum);

% Create contour plot
subplot(1,2,1);
contourf(Departure_dn, Arrival_dn, delta_v, 50, 'LineStyle', 'none');
colorbar;
clim([min(delta_v) max(delta_v)]);
xlabel('Departure Date');
ylabel('Arrival Date');
title('Earth-'+ planetName + ' Porkchop Plot (\DeltaV, km/s)');
grid on;

% Add contour lines
hold on;
[C2, h2] = contour(Departure_dn, Arrival_dn, delta_v, 15, 'k', 'LineWidth', 0.5);
clabel(C2, h2, 'FontSize', 8, 'Color', 'white');

% Format dates on axes
datetick('x', 'mm/dd/yy', 'keeplimits');
datetick('y', 'mm/dd/yy', 'keeplimits');

% Find and mark optimal transfer
valid_dv = delta_v;
valid_dv(isnan(valid_dv)) = inf;
[min_dv, min_idx] = min(valid_dv(:));

if isfinite(min_dv)
    [opt_arr_idx, opt_dep_idx] = ind2sub(size(delta_v), min_idx);
    opt_dep_datenum = departure_datenum(opt_dep_idx);
    opt_arr_datenum = arrival_datenum(opt_arr_idx);
    
    plot(opt_dep_datenum, opt_arr_datenum, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'red');
    text(opt_dep_datenum, opt_arr_datenum, sprintf('  Optimal: %.2f km/s', min_dv), ...
         'Color', 'red', 'FontWeight', 'bold');
end

%% Time of Flight Contour
subplot(1,2,2);

% Calculate time of flight in days
tof_days = Arrival - Departure;

% Create ToF contour
contourf(Departure_dn, Arrival_dn, tof_days, 20, 'LineStyle', 'none');
colorbar;
xlabel('Departure Date');
ylabel('Arrival Date');
title('Time of Flight (days)');
grid on;

% Add contour lines
hold on;
[C_tof2, h_tof2] = contour(Departure_dn, Arrival_dn, tof_days, 10, 'k', 'LineWidth', 0.5);
clabel(C_tof2, h_tof2, 'FontSize', 8, 'Color', 'white');

% Mark optimal transfer
if isfinite(min_dv)
    plot(opt_dep_datenum, opt_arr_datenum, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'red');
end

datetick('x', 'mm/dd/yy', 'keeplimits');
datetick('y', 'mm/dd/yy', 'keeplimits');

%% Display optimal transfer information
if isfinite(min_dv)
    fprintf('\n=== Optimal Transfer Found ===\n');
    opt_dep_date = datetime(opt_dep_datenum, 'ConvertFrom', 'datenum');
    opt_arr_date = datetime(opt_arr_datenum, 'ConvertFrom', 'datenum');
    fprintf('Departure Date: %s\n', datestr(opt_dep_date, 'mmm dd, yyyy, HH:MM:SS'));
    fprintf('Arrival Date: %s\n', datestr(opt_arr_date, 'mmm dd, yyyy, HH:MM:SS'));
    fprintf('Time of Flight: %.1f days\n', tof_days(opt_arr_idx, opt_dep_idx));
    fprintf('Total ΔV: %.3f km/s\n', min_dv);
else
    fprintf('\nNo valid transfers found in the specified window.\n');
end


%% Create Planet as Spheres

% Create the sphere
[XS,YS,ZS] = sphere(50);

% Scale the sphere of Earth
XS_Earth = XS * R_E;
YS_Earth = YS * R_E;
ZS_Earth = ZS * -R_E;


% Scale the sphere of Planet
XS_Planet = XS * R_planet;
YS_Planet = YS * R_planet;
ZS_Planet = ZS * -R_planet;


% Scale the sphere of Sun
XS_Sun = XS * R_Sun;
YS_Sun = YS * R_Sun;
ZS_Sun = ZS * -R_Sun;


%% Plot Optimal Transfer Orbit (IMPROVED VERSION)
if isfinite(min_dv)
    fprintf('Plotting optimal transfer orbit...\n');
    
    % Get optimal dates
    opt_dep_mjd2000 = departure_dates(opt_dep_idx);
    opt_arr_mjd2000 = arrival_dates(opt_arr_idx);
    opt_tof_days = opt_arr_mjd2000 - opt_dep_mjd2000;
    opt_tof_seconds = opt_tof_days * 24 * 3600;
    
    % Get planetary positions for optimal transfer
    [kep_earth_opt, ~] = uplanet(opt_dep_mjd2000, 3);
    kep_earth_deg_opt = [kep_earth_opt(1),kep_earth_opt(2),rad2deg(kep_earth_opt(3)),rad2deg(kep_earth_opt(4)), ...
        rad2deg(kep_earth_opt(5)),rad2deg(kep_earth_opt(6))];
    
    [r1_opt, v1_opt] = kep2car(kep_earth_deg_opt, mu_Sun);

    y_earth = [r1_opt; v1_opt];
    T_earth = 2*pi*sqrt(kep_earth_deg_opt(1)^3/mu_Sun); % Orbital period [s]
    tspan_earth = linspace( 0, T_earth, 10000); % Reduced to 10 periods for clarity

    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

    % Perform the integration
    [ Tearth, Yearth ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_earth, y_earth, options );

    % Calculate ground track
    r_earth = Yearth(:,1:3);
    %gtrack_data = zeros(length(T),4);

    % In your main script, before the ground track loop:
    %fprintf('First position vector: [%f, %f, %f]\n', r_2(1,1), r_2(1,2), r_2(1,3));
    %fprintf('Time range: %f to %f seconds\n', tspan_2(1), tspan_2(end));


    
    [kep_mars_opt, ~] = uplanet(opt_arr_mjd2000, N_planet);
    kep_mars_deg_opt = [kep_mars_opt(1),kep_mars_opt(2),rad2deg(kep_mars_opt(3)),rad2deg(kep_mars_opt(4)), ...
        rad2deg(kep_mars_opt(5)),rad2deg(kep_mars_opt(6))];
    [r2_opt, v2_opt] = kep2car(kep_mars_opt, mu_Sun);

    y_planet = [r2_opt; v2_opt];
    T_planet = 2*pi*sqrt(kep_mars_deg_opt(1)^3/mu_Sun); % Orbital period [s]
    tspan_planet = linspace( 0, T_planet, 10000);

    % Set options for the ODE solver
    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

    % Perform the integration
    [ Tplanet, Yplanet ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_planet, y_planet, options );

    % Calculate ground track
    r_planet = Yplanet(:,1:3);
    
    % Ensure vectors are column vectors
    r1_opt = r1_opt(:);
    v1_opt = v1_opt(:);
    r2_opt = r2_opt(:);
    v2_opt = v2_opt(:);
    
    % Solve Lambert for optimal transfer
    [A_opt, P_opt, E_opt, ERROR_opt, VI_opt, VF_opt, ~, THETA_opt] = ...
        lambertMR(r1_opt, r2_opt, opt_tof_seconds, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);

    VI_opt = VI_opt(:);
    VF_opt = VF_opt(:);
    
    if ERROR_opt == 0
        % Ensure velocity vectors are column vectors
        VI_opt = VI_opt(:);
        VF_opt = VF_opt(:);
        
        % Create figure for optimal transfer
        figure(3)
        
        % 3D view

        subplot(2,3,[1,2,4,5]);
        hold on; grid on; axis equal;
        
        % Plot Sun with proper scaling
        surf(XS_Sun, YS_Sun, ZS_Sun, 'FaceColor', 'texturemap', 'CData', SunImage, 'EdgeColor', 'none','DisplayName','Sun');
        
        % Generate planetary orbits using analytical approach (more reliable)
        fprintf('Generating planetary orbits...\n');
        
        
        % Plot planetary orbits
        plot3(r_earth(:,1), r_earth(:,2), r_earth(:,3), 'b-', 'LineWidth', 0.5, 'DisplayName', 'Earth Orbit');
        plot3(r_planet(:,1), r_planet(:,2), r_planet(:,3), 'r-', 'LineWidth', 0.5, 'DisplayName', planetName + ' Orbit');
        
        % Generate and plot transfer orbit using propagation
        fprintf('Propagating transfer orbit...\n');
        y0_transfer = [r1_opt; VI_opt];
        
        try
            [~, Y_transfer] = ode113(@(t,y) ode_2bp(t, y, mu_Sun), ...
                                    linspace(0, opt_tof_seconds, 1000), y0_transfer);
            plot3(Y_transfer(:,1), Y_transfer(:,2), Y_transfer(:,3), 'g-', 'LineWidth', 2.5, ...
                  'DisplayName', 'Transfer Orbit');
        catch
            fprintf('ODE propagation failed, using analytical transfer orbit.\n');
            % Fallback: plot direct connection
            plot3([r1_opt(1), r2_opt(1)], [r1_opt(2), r2_opt(2)], [r1_opt(3), r2_opt(3)], ...
                  'g-', 'LineWidth', 2.5, 'DisplayName', 'Transfer Orbit');
        end
        
        % Mark key positions
        surf(XS_Earth + r1_opt(1), YS_Earth + r1_opt(2), ZS_Earth + r1_opt(3), 'FaceColor', 'texturemap', 'CData', EarthImage, 'EdgeColor', 'none','DisplayName','Earth');
        surf(XS_Planet + r2_opt(1), YS_Planet + r2_opt(2), ZS_Planet + r2_opt(3), 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none','DisplayName',planetName);
        
        % Plot velocity vectors (scaled for visibility)
        scale = 1e3; % Adjust this scaling factor as needed
        if length(VI_opt) == 3
            quiver3(r1_opt(1), r1_opt(2), r1_opt(3), ...
                    VI_opt(1)*scale, VI_opt(2)*scale, VI_opt(3)*scale, ...
                    'b-', 'LineWidth', 0.5, 'MaxHeadSize', 1, 'DisplayName', 'Departure Velocity');
        end
        
        if length(VF_opt) == 3
            quiver3(r2_opt(1), r2_opt(2), r2_opt(3), ...
                    VF_opt(1)*scale, VF_opt(2)*scale, VF_opt(3)*scale, ...
                    'r-', 'LineWidth', 0.5, 'MaxHeadSize', 1, 'DisplayName', 'Arrival Velocity');
        end
        
        xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
        title(sprintf('Optimal Earth-'+ planetName + ' Transfer Orbit\nDeparture: %s, Arrival: %s\n\\DeltaV: %.3f km/s, ToF: %.1f days', ...
              datestr(datetime(opt_dep_mjd2000 + 730486.5, 'ConvertFrom', 'juliandate'), 'mmm dd, yyyy'), ...
              datestr(datetime(opt_arr_mjd2000 + 730486.5, 'ConvertFrom', 'juliandate'), 'mmm dd, yyyy'), ...
              min_dv, opt_tof_days));
        legend('Location', 'best');
        view(45, 30);
        
        % XY plane projection
        subplot(2,3,3);
        hold on; grid on; axis equal;
        
        % Plot orbits in XY plane
        plot(r_earth(:,1), r_earth(:,2), 'b-', 'LineWidth', 1);
        plot(r_planet(:,1), r_planet(:,2), 'r-', 'LineWidth', 1);
        if exist('Y_transfer', 'var')
            plot(Y_transfer(:,1), Y_transfer(:,2), 'g-', 'LineWidth', 2);
        else
            plot([r1_opt(1), r2_opt(1)], [r1_opt(2), r2_opt(2)], 'g-', 'LineWidth', 2);
        end
        
        % Plot Sun
        plot(0, 0, 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
        
        % Mark positions
        plot(r1_opt(1), r1_opt(2), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'blue');
        plot(r2_opt(1), r2_opt(2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
        
        xlabel('X (km)'); ylabel('Y (km)');
        title('XY Plane Projection');
        
        % XZ plane projection
        subplot(2,3,6);
        hold on; grid on; axis equal;
        
        % Plot orbits in XZ plane
        plot(r_earth(:,1), r_earth(:,3), 'b-', 'LineWidth', 1);
        plot(r_planet(:,1), r_planet(:,3), 'r-', 'LineWidth', 1);
        if exist('Y_transfer', 'var')
            plot(Y_transfer(:,1), Y_transfer(:,3), 'g-', 'LineWidth', 2);
        else
            plot([r1_opt(1), r2_opt(1)], [r1_opt(3), r2_opt(3)], 'g-', 'LineWidth', 2);
        end
        
        % Plot Sun
        plot(0, 0, 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
        
        % Mark positions
        plot(r1_opt(1), r1_opt(3), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'blue');
        plot(r2_opt(1), r2_opt(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
        
        xlabel('X (km)'); ylabel('Z (km)');
        title('XZ Plane Projection');
        
        % Display detailed transfer parameters
        fprintf('\n=== Detailed Transfer Parameters ===\n');
        fprintf('Earth position at departure: [%.2f, %.2f, %.2f] km\n', r1_opt);
        fprintf(planetName + ' position at arrival: [%.2f, %.2f, %.2f] km\n', r2_opt);
        fprintf('Earth velocity: [%.4f, %.4f, %.4f] km/s\n', v1_opt);
        fprintf( planetName +' velocity: [%.4f, %.4f, %.4f] km/s\n', v2_opt);
        fprintf('Transfer orbit velocity at Earth: [%.4f, %.4f, %.4f] km/s\n', VI_opt);
        fprintf('Transfer orbit velocity at '+ planetName + ': [%.4f, %.4f, %.4f] km/s\n', VF_opt);
        fprintf('Semi-major axis: %.6f AU\n', A_opt / AU);
        fprintf('Eccentricity: %.6f\n', E_opt);
        fprintf('Transfer angle: %.2f°\n', rad2deg(THETA_opt));
        fprintf('Departure ΔV: %.3f km/s\n', norm(VI_opt - v1_opt));
        fprintf('Arrival ΔV: %.3f km/s\n', norm(VF_opt - v2_opt));
        
    else
        fprintf('Error computing optimal transfer orbit: %d\n', ERROR_opt);
    end
    
else
    fprintf('No valid optimal transfer found to plot.\n');
end

fprintf('\nAnalysis completed successfully!\n');