%% Assignment1.m - Corrected version with all fixes applied
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
n_points = 50;  % Increase to 100+ after verification
% Departure window from Earth (2030-2060)
t_dep_e = date2mjd2000([2031, 1, 1, 0, 0, 0]);
t_dep_l = date2mjd2000([2032, 1, 1, 0, 0, 0]);
departure_dates = linspace(t_dep_e, t_dep_l, n_points);
% Gravity assist at Mars
t_ga_e = date2mjd2000([2033, 1, 1, 0, 0, 0]);
t_ga_l = date2mjd2000([2034, 1, 1, 0, 0, 0]);
ga_dates = linspace(t_ga_e, t_ga_l, n_points);
% Arrival at asteroid
t_arr_e = date2mjd2000([2035, 1, 1, 0, 0, 0]);
t_arr_l = date2mjd2000([2036, 1, 1, 0, 0, 0]);
arrival_dates = linspace(t_arr_e, t_arr_l, n_points);
% Initialize matrices
delta_v_total = NaN(n_points, n_points, n_points);
flyby_feasible = false(n_points, n_points, n_points);
min_periapsis = NaN(n_points, n_points, n_points);
fprintf('Computing porkchop plot...\n');
progress_interval = max(1, floor(n_points^3/100));
counter = 0;
tof1days = zeros(n_points,n_points,n_points);
tof2days = zeros(n_points,n_points,n_points);
dv_flyby = zeros(n_points,n_points,n_points);
dv_departure = zeros(n_points,n_points,n_points);
dv_arrival = zeros(n_points,n_points,n_points);
%% Main Computation Loop - CORRECTED
for i = 1:n_points      
    for j = 1:n_points  
        for k = 1:n_points  
            counter = counter + 1;
            if mod(counter, progress_interval) == 0
                fprintf('Progress: %.1f%%\n', (counter/n_points^3)*100);
            end
            t_dep = departure_dates(i);
            t_ga = ga_dates(j);
            t_arr = arrival_dates(k);
            tof1days(i,j,k) = t_ga - t_dep;
            tof2days(i,j,k) = t_arr - t_ga;
            if ~(t_dep < t_ga && t_ga < t_arr)
                continue;
            end
            tof_depga = (t_ga - t_dep) * 86400;  
            tof_gaarr = (t_arr - t_ga) * 86400;  
            try
                % Get positions/velocities - ALL COLUMN VECTORS
                kep_earth = uplanet(t_dep, 3);
                [r_earth, v_earth] = kep2car(kep_earth, mu_Sun);
                r_earth = r_earth(:); v_earth = v_earth(:);
                
                kep_mars = uplanet(t_ga, 4);    
                [r_mars, v_mars] = kep2car(kep_mars, mu_Sun);
                r_mars = r_mars(:); v_mars = v_mars(:);
                
                kep_aster = ephAsteroids(t_arr, 257323);
                [r_aster, v_aster] = kep2car(kep_aster, mu_Sun);
                r_aster = r_aster(:); v_aster = v_aster(:);  % FIX: Added
                
                % Lambert arcs
                [~, ~, ~, error1, VI1, VF1, ~, ~] = lambertMR(r_earth, r_mars, tof_depga, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);
                [~, ~, ~, error2, VI2, VF2, ~, ~] = lambertMR(r_mars, r_aster, tof_gaarr, mu_Sun, orbitType, Nrev, Ncase, optionsLMR);
                
                VI1 = VI1(:); VF1 = VF1(:); VI2 = VI2(:); VF2 = VF2(:);
                
                if error1 ~= 0 || error2 ~= 0
                    continue;
                end
                
                % FIXED: v_inf definitions (relative to planet)
                v_inf_minus = VF1 - v_mars;     % Incoming v_inf at Mars
                v_inf_plus  = VI2 - v_mars;     % Outgoing v_inf from Mars
                
                v_in_norm = norm(v_inf_minus);
                v_out_norm = norm(v_inf_plus);
                
                if v_in_norm == 0 || v_out_norm == 0
                    continue;
                end
                
                % Turning angle
                cos_delta = dot(v_inf_minus, v_inf_plus) / (v_in_norm * v_out_norm);
                cos_delta = max(-1, min(1, cos_delta));
                delta = acos(cos_delta);
                
                if ~(delta > 0 && delta < pi)
                    continue;
                end
                
                % FIXED fzero: Bounded solver + better checks
                safety_margin = 200;
                r_p_min = R_Mars + safety_margin;
                r_p_max = 5e4;  % Reasonable max altitude
                
                fun = @(r_p) asin(1/(1 + r_p*v_out_norm^2/mu_Mars)) + ...
                             asin(1/(1 + r_p*v_in_norm^2/mu_Mars)) - delta;
                
                try
                    options = optimset('TolX', 1e-12, 'TolFun', 1e-16, 'MaxIter', 200, 'Display', 'off');
                    r_p = fzero(fun, [r_p_min, r_p_max], options);  % Bounded!
                    
                    if isfinite(r_p) && r_p >= r_p_min && fun(r_p) < 1e-10  % Verify root
                        min_periapsis(i,j,k) = r_p;
                        flyby_feasible(i,j,k) = true;
                        
                        % FIXED dv_flyby: Correct powered flyby delta-v
                        % dv_flyby(i,j,k) = abs(v_out_norm - v_in_norm);

                        % Powered flyby ΔV computed at periapsis (NOT at v_inf)
                        v_p_in  = sqrt(v_in_norm^2  + 2*mu_Mars/r_p);
                        v_p_out = sqrt(v_out_norm^2 + 2*mu_Mars/r_p);
                        dv_flyby(i,j,k) = abs(v_p_out - v_p_in);

                        
                        dv_departure(i,j,k) = norm(VI1 - v_earth);
                        dv_arrival(i,j,k) = norm(VF2 - v_aster);  % Now safe dims
                        
                        % Total - now finite!
                        delta_v_total(i,j,k) = dv_departure(i,j,k) + dv_flyby(i,j,k) + dv_arrival(i,j,k);
                        
                    end
                catch
                    % Silent fail
                end
                
            catch ME
                % Silent fail on errors
            end
        end
    end
end
fprintf('Computation complete!\n');
%% Find Optimal Trajectory - CORRECTED
feasible_mask = flyby_feasible & isfinite(delta_v_total);
feasible_indices = find(feasible_mask);
if ~isempty(feasible_indices)
    [min_dv, min_idx] = min(delta_v_total(feasible_indices));
    [i_opt, j_opt, k_opt] = ind2sub(size(delta_v_total), feasible_indices(min_idx));
    
    t_dep_opt = departure_dates(i_opt);
    t_ga_opt = ga_dates(j_opt);
    t_arr_opt = arrival_dates(k_opt);
    
    % Dates
    dep_date_vec = mjd20002date(t_dep_opt);
    ga_date_vec = mjd20002date(t_ga_opt);
    arr_date_vec = mjd20002date(t_arr_opt);
    
    dep_datetime = datetime(dep_date_vec(1), dep_date_vec(2), dep_date_vec(3), ...
                           dep_date_vec(4), dep_date_vec(5), floor(dep_date_vec(6)));
    ga_datetime = datetime(ga_date_vec(1), ga_date_vec(2), ga_date_vec(3), ...
                          ga_date_vec(4), ga_date_vec(5), floor(ga_date_vec(6)));
    arr_datetime = datetime(arr_date_vec(1), arr_date_vec(2), arr_date_vec(3), ...
                           arr_date_vec(4), arr_date_vec(5), floor(ga_date_vec(6)));
    
    fprintf('\n=== OPTIMAL TRAJECTORY FOUND ===\n');
    fprintf('Departure: %s (MJD2000 = %.2f)\n', string(dep_datetime), t_dep_opt);
    fprintf('Gravity Assist: %s (MJD2000 = %.2f)\n', string(ga_datetime), t_ga_opt);
    fprintf('Arrival: %s (MJD2000 = %.2f)\n', string(arr_datetime), t_arr_opt);
    fprintf('Total ΔV: %.4f km/s\n', min_dv);
    fprintf('ΔV Departure: %.4f km/s\n', dv_departure(i_opt,j_opt,k_opt));
    fprintf('ΔV Powered Gravity Assist: %.4f km/s\n', dv_flyby(i_opt,j_opt,k_opt));
    fprintf('ΔV Arrival: %.4f km/s\n', dv_arrival(i_opt,j_opt,k_opt));
    fprintf('Time of Flight: %.2f days\n', t_arr_opt - t_dep_opt);
    fprintf('Flyby periapsis: %.2f km (%.2f R_Mars)\n', ...
            min_periapsis(i_opt,j_opt,k_opt), ...
            min_periapsis(i_opt,j_opt,k_opt)/R_Mars);
    
    fprintf('\n=== TRAJECTORY DETAILS ===\n');
    fprintf('Earth-Mars transfer time: %.2f days\n', t_ga_opt - t_dep_opt);
    fprintf('Mars-Asteroid transfer time: %.2f days\n', t_arr_opt - t_ga_opt);
    fprintf('Flyby altitude: %.2f km\n', min_periapsis(i_opt,j_opt,k_opt) - R_Mars);        
    
    % FIXED: Recalculate trajectory data with correct vectors
    tof_depga = (t_ga_opt - t_dep_opt) * 86400;
    tof_gaarr = (t_arr_opt - t_ga_opt) * 86400;
    
    kep_earth = uplanet(t_dep_opt, 3); [r_earth, v_earth] = kep2car(kep_earth, mu_Sun); r_earth=r_earth(:); v_earth=v_earth(:);
    kep_mars = uplanet(t_ga_opt, 4); [r_mars, v_mars] = kep2car(kep_mars, mu_Sun); r_mars=r_mars(:); v_mars=v_mars(:);
    kep_aster = ephAsteroids(t_arr_opt, 257323); [r_aster, v_aster] = kep2car(kep_aster, mu_Sun); r_aster=r_aster(:); v_aster=v_aster(:);
    
    [~, ~, ~, ~, VI1, VF1, ~, ~] = lambertMR(r_earth, r_mars, tof_depga, mu_Sun, 0, 0, 0, 0);
    [~, ~, ~, ~, VI2, VF2, ~, ~] = lambertMR(r_mars, r_aster, tof_gaarr, mu_Sun, 0, 0, 0, 0);
    
    VI1=VI1(:); VF1=VF1(:); VI2=VI2(:); VF2=VF2(:);
    
    % FIXED v_inf - no transpose
    v_inf_minus = VF1 - v_mars;
    v_inf_plus = VI2 - v_mars;
    
    v_in_norm = norm(v_inf_minus);
    v_out_norm = norm(v_inf_plus);
    delta = acos( dot(v_inf_minus, v_inf_plus) / (v_in_norm * v_out_norm) );
    
    % Store corrected trajectory
    opt_trajectory.t_dep = t_dep_opt; opt_trajectory.t_ga = t_ga_opt; opt_trajectory.t_arr = t_arr_opt;
    opt_trajectory.r_earth = r_earth; opt_trajectory.r_mars = r_mars; opt_trajectory.r_aster = r_aster;
    opt_trajectory.VI1 = VI1; opt_trajectory.VF1 = VF1; opt_trajectory.VI2 = VI2; opt_trajectory.VF2 = VF2;
    opt_trajectory.v_earth = v_earth; opt_trajectory.v_mars = v_mars; opt_trajectory.v_aster = v_aster;
    
    opt_flyby.r_p = min_periapsis(i_opt,j_opt,k_opt);
    opt_flyby.v_inf_minus = v_inf_minus; opt_flyby.v_inf_plus = v_inf_plus;
    opt_flyby.v_in_norm = v_in_norm; opt_flyby.v_out_norm = v_out_norm;
    opt_flyby.delta = delta;
    
else
    fprintf('\nNo feasible trajectories found!\n');
    return;
end
%% Plotting Section - CORRECTED
fprintf('\nGenerating plots...\n');
figure('Name', 'Porkchop Plot');
subplot(1, 3, 1);
dv_slice_dep_ga = squeeze(dv_departure(:,:,k_opt));
create_porkchop_plot(departure_dates, ga_dates, dv_slice_dep_ga, i_opt, j_opt, tof1days(:,:,k_opt), dv_departure(i_opt,j_opt,k_opt));
title('Departure vs Gravity Assist (Fixed Arrival)');

subplot(1, 3, 2);
dv_slice_ga_arr = squeeze(dv_arrival(i_opt,:,:));
create_porkchop_plot(ga_dates, arrival_dates, dv_slice_ga_arr, j_opt, k_opt, tof2days(i_opt,:,:), dv_arrival(i_opt,j_opt,k_opt));
title('Gravity Assist vs Arrival (Fixed Departure)');

subplot(1, 3, 3);
dv_slice_dep_arr = squeeze(delta_v_total(:,j_opt,:));
create_porkchop_plot(departure_dates, arrival_dates, dv_slice_dep_arr, i_opt, k_opt, (tof1days(:,j_opt,:) + tof2days(:,j_opt,:)), delta_v_total(i_opt,j_opt,k_opt));
title('Departure vs Arrival (Fixed GA)');

%% Comprehensive Plots - BIG 1x2 LAYOUT
figure('Name', 'Optimized Trajectory Analysis', 'Units','normalized', 'Position',[0.05 0.08 0.90 0.84]);

tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');  % recommended over subplot [page:63]

ax3d = nexttile(1);
plot_3d_trajectory_comprehensive(opt_trajectory, opt_flyby, SunImage);
title(ax3d, '3D View: Complete Transfer Orbit');

ax2d = nexttile(2);
plot_2d_projection(opt_trajectory);
title(ax2d, '2D Projection: Ecliptic Plane View');

% sgtitle('Optimized Trajectory: Earth \rightarrow Mars \rightarrow Asteroid (257323)', ...
%         'FontSize', 14, 'FontWeight', 'bold');
% sgt.Units = 'normalized';
% sgt.Position(2) = sgt.Position(2) - 0.09;   % move up (tune 0.02–0.06)
%% CORRECTED Helper Functions
function plot_3d_trajectory_comprehensive(traj, ~, img1)
    mu_Sun = astroConstants(4);
    R_Sun  = astroConstants(3);

    hold on; grid on; box on; axis equal; view(45, 25);

    % Sun
    [Xs,Ys,Zs] = sphere(60);
    if ~isempty(img1)
        surf(25*R_Sun*Xs, 25*R_Sun*Ys, 25*R_Sun*Zs, ...
            'FaceColor','texturemap','CData',img1,'EdgeColor','none','FaceAlpha',0.85, ...
            'DisplayName','Sun');
    else
        surf(25*R_Sun*Xs, 25*R_Sun*Ys, 25*R_Sun*Zs, ...
            'FaceColor',[1 1 0],'EdgeColor','none','FaceAlpha',0.5, ...
            'DisplayName','Sun');
    end

    % 1) Earth orbit
    plot_planet_orbit(traj.r_earth, traj.v_earth, mu_Sun, 'b', 'Earth orbit');

    % 2) Mars orbit
    plot_planet_orbit(traj.r_mars, traj.v_mars, mu_Sun, 'r', 'Mars orbit');

    % 3) Asteroid orbit
    plot_planet_orbit(traj.r_aster, traj.v_aster, mu_Sun, [0.3 1 0.3], 'Asteroid orbit');

    % Epoch markers
    plot3(traj.r_earth(1), traj.r_earth(2), traj.r_earth(3), 'bo', ...
        'MarkerFaceColor','b','DisplayName','Earth @ dep');
    plot3(traj.r_mars(1), traj.r_mars(2), traj.r_mars(3), 'ro', ...
        'MarkerFaceColor','r','DisplayName','Mars @ GA');
    plot3(traj.r_aster(1), traj.r_aster(2), traj.r_aster(3), 'go', ...
        'MarkerFaceColor','g','DisplayName','Asteroid @ arr');

    % 4) Lambert arc 1
    plot_transfer_arc(traj.r_earth, traj.VI1, traj.t_ga - traj.t_dep, mu_Sun, 'c', 'Transfer 1 (Earth→Mars)');

    % 5) Lambert arc 2 (must use VI2, not VF2)
    plot_transfer_arc(traj.r_mars, traj.VI2, traj.t_arr - traj.t_ga, mu_Sun, 'm', 'Transfer 2 (Mars→Asteroid)');

    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

    max_dist = 1.10 * max([norm(traj.r_earth), norm(traj.r_mars), norm(traj.r_aster)]);
    xlim([-max_dist max_dist]);
    ylim([-max_dist max_dist]);

    % Don't crush Z: keep symmetric range (still smaller than XY because near-coplanar)
    zlim([-0.25*max_dist 0.25*max_dist]);

    lgd = legend('Location','bestoutside');
    lgd.AutoUpdate = 'off';

    camlight('headlight'); lighting gouraud;
    hold off;
end



function plot_2d_projection(traj)
    mu_Sun = astroConstants(4);
    hold on; grid on; box on; axis equal;
    plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'y');  % Sun
    
    % Orbits 2D
    plot_planet_orbit_2d(traj.r_earth, traj.v_earth, mu_Sun, 'b');
    plot_planet_orbit_2d(traj.r_mars, traj.v_mars, mu_Sun, 'r');
    plot_planet_orbit_2d(traj.r_aster, traj.v_aster, mu_Sun, [0.3 1 0.3]);
    
    % FIXED arcs 2D
    plot_transfer_arc_2d(traj.r_earth, traj.VI1, traj.t_ga - traj.t_dep, mu_Sun, 'b');
    plot_transfer_arc_2d(traj.r_mars, traj.VI2, traj.t_arr - traj.t_ga, mu_Sun, 'm');
    
    % Positions
    plot(traj.r_earth(1), traj.r_earth(2), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    plot(traj.r_mars(1), traj.r_mars(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(traj.r_aster(1), traj.r_aster(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    
    legend('Location', 'best'); xlabel('X (km)'); ylabel('Y (km)');
    max_dist = 1.1 * max([norm(traj.r_earth(1:2)), norm(traj.r_mars(1:2)), norm(traj.r_aster(1:2))]);
    xlim([-max_dist max_dist]); ylim([-max_dist max_dist]); hold off;
end

function plot_transfer_arc(r0, v0, tof_days, mu, color, name)
    tof_sec = tof_days * 86400;
    tspan = linspace(0, tof_sec, 1000);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);
    y0 = [r0; v0];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options);
    plot3(Y(:,1), Y(:,2), Y(:,3), color, 'LineWidth', 2, 'DisplayName', name);
end

function plot_transfer_arc_2d(r0, v0, tof_days, mu, color)
    tof_sec = tof_days * 86400;
    tspan = linspace(0, tof_sec, 1000);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);

    r0 = r0(:); v0 = v0(:);
    y0 = [r0; v0];

    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options);

    plot(Y(:,1), Y(:,2), 'Color', color, 'LineWidth', 2);
end


function plot_planet_orbit(r0, v0, mu, color, name)
    r0 = r0(:); v0 = v0(:);

    kep = car2kep(r0, v0, mu);
    a = kep(1);
    T = 2*pi*sqrt(a^3/mu);

    tspan = linspace(0, T, 5000);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);

    y0 = [r0; v0];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options);

    % IMPORTANT: use 'Color', not passing RGB as a positional linespec
    plot3(Y(:,1), Y(:,2), Y(:,3), 'Color', color, 'LineWidth', 0.8, 'DisplayName', name);
end


function plot_planet_orbit_2d(r0, v0, mu, color)
    plot_planet_orbit(r0, v0, mu, color, '');
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
