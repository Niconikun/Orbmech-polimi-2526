%% Exercise 3 – Mars Express (porkchop + trajectory)
clear; clc; close all;
addpath("Functions 2bp\")
addpath("time\")

%% ==== CONSTANTS & HELPERS ====
deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;

% Mission windows (from slides)
% Departure:  2003 April 1 – 2003 August 1
% Arrival:    2003 September 1 – 2004 March 1
dep_start = datenum(2024,6,1,0,0,0);
dep_end   = datenum(2026,11,1,0,0,0);
arr_start = datenum(2024,12,1,0,0,0);
arr_end   = datenum(2027,6,1,0,0,0);

max_launcher = 7; %km/s

% Grid sizes (increase for smoother plots)
N_dep = 700;   % number of departure samples
N_arr = 700;   % number of arrival samples

dep_vec = linspace(dep_start, dep_end, N_dep);  % MATLAB datenums
arr_vec = linspace(arr_start, arr_end, N_arr);

% Output arrays: rows = arrival, cols = departure
DVtot   = NaN(N_arr, N_dep);    % total Delta-v [km/s]
TOFdays = NaN(N_arr, N_dep);    % time of flight [days]

% For Lambert
orbitType  = 0;   % prograde
Nrev       = 0;   % zero-revolution
Ncase      = 0;
optionsLMR = 0;   % no display

%% ==== BUILD PORKCHOP (DOUBLE LOOP) ====
fprintf('Building porkchop grid (%d x %d)...\n', N_arr, N_dep);

muSun = [];  % we will get it from uplanet the first time

for iA = 1:N_arr      % arrival index (rows)
    for jD = 1:N_dep  % departure index (cols)

        dep_num = dep_vec(jD);   % MATLAB datenum
        arr_num = arr_vec(iA);

        tof_days = arr_num - dep_num;
        if tof_days <= 0
            continue;    % invalid combination
        end

        % Convert to MJD2000 for uplanet
        [yd, md, dd, hd, mind, sd] = datevec(dep_num);
        [ya, ma, da, ha, mina, sa] = datevec(arr_num);

        mjd_dep = date2mjd2000_local(yd, md, dd, hd, mind, sd);
        mjd_arr = date2mjd2000_local(ya, ma, da, ha, mina, sa);

        % Earth & Mars heliocentric states at dep & arr
        [kepE, ksun_dep] = uplanet(mjd_dep, 1);   % Earth
        [kepM, ~]        = uplanet(mjd_arr, 2);   % Mars

        if isempty(muSun)
            muSun = ksun_dep;   % [km^3/s^2]
        end

        [rE, vE] = kep2car_local(kepE, muSun);
        [rM, vM] = kep2car_local(kepM, muSun);

        tof_sec = day2sec(tof_days);

        % Lambert (Earth -> Mars)
        [A, P, Ecc, ERROR, v1_t, v2_t, Tpar, THETA] = ...
            lambertMR(rE.', rM.', tof_sec, muSun, orbitType, Nrev, Ncase, optionsLMR);

        if ERROR ~= 0
            % no solution / convergence, leave NaN
            continue;
        end

        dv1 = norm(v1_t - vE.');   % departure hyperbolic injection (heliocentric)
        dv2 = norm(vM.' - v2_t);   % arrival insertion (heliocentric)
        if dv1 <= max_launcher
            dv_tot = dv1 + dv2;  % Calculate total Delta-v
        else
            dv_tot = NaN;  % Set total Delta-v to NaN if the launcher limit is exceeded
        end

        DVtot(iA,jD)   = dv_tot;
        TOFdays(iA,jD) = tof_days;
    end
end

%% ==== FIND MINIMUM Δv IN GRID ====
[dv_min, idxMin] = min(DVtot(:), [], 'omitnan');
[iA_min, jD_min] = ind2sub(size(DVtot), idxMin);

dep_opt_num = dep_vec(jD_min);
arr_opt_num = arr_vec(iA_min);
tof_opt_days = TOFdays(iA_min, jD_min);

fprintf('\nGRID MINIMUM:\n');
fprintf('  Δv_min   = %.4f km/s\n', dv_min);
fprintf('  Departure = %s\n', datestr(dep_opt_num, 31));
fprintf('  Arrival   = %s\n', datestr(arr_opt_num, 31));
fprintf('  ToF       = %.2f days\n\n', tof_opt_days);

%% ==== PORKCHOP PLOT ====
[DEPgrid, ARRgrid] = meshgrid(dep_vec, arr_vec);

figure;
% Δv contours (porkchop)
% instead of contourf(...,30,...)
contourf(DEPgrid, ARRgrid, DVtot, linspace(5,10,30), 'LineStyle','none');
colormap jet;
c = colorbar;
caxis([5 10]);                 % *<- key line*
c.Ticks = 5:0.5:10;
c.Label.String = '\Delta v [km/s]';

hold on; grid on;
colorbar;
colormap jet;
xlabel('Departure date');
ylabel('Arrival date');
title('\DeltaV_{tot} porkchop – Mars Express');

% Format axes as calendar dates
datetick('x','dd-mmm-03','keeplimits');
datetick('y','dd-mmm-03','keeplimits');

% Time-of-flight contour lines (days)
ToF_levels = 100:20:300; % choose TOF isolines [days]
[C2,h2] = contour(DEPgrid, ARRgrid, TOFdays, ToF_levels, 'k--', 'LineWidth', 0.8);
clabel(C2,h2,'Color','k','FontSize',8);

% Mark grid minimum
plot(dep_opt_num, arr_opt_num, 'wo', 'MarkerFaceColor','k', 'MarkerSize',7);
text(dep_opt_num, arr_opt_num, sprintf('  min \\DeltaV = %.2f km/s', dv_min), ...
    'Color','w','FontSize',8,'FontWeight','bold');

hold off;

%% ==== TRANSFER TRAJECTORY FOR GRID MINIMUM ====
fprintf('Computing trajectory for grid minimum...\n');

% States at optimal departure/arrival (ephemerides)
[yd, md, dd, hd, mind, sd] = datevec(dep_opt_num);
[ya, ma, da, ha, mina, sa] = datevec(arr_opt_num);

mjd_dep_opt = date2mjd2000_local(yd, md, dd, hd, mind, sd);
mjd_arr_opt = date2mjd2000_local(ya, ma, da, ha, mina, sa);

[kepE_dep, muSun] = uplanet(mjd_dep_opt, 1);
[kepM_arr, ~]     = uplanet(mjd_arr_opt, 2);

[rE_dep, vE_dep] = kep2car_local(kepE_dep, muSun);
[rM_arr, vM_arr] = kep2car_local(kepM_arr, muSun);

tof_opt_sec = day2sec(tof_opt_days);

[~, ~, ~, ERRORopt, v1_t_opt, v2_t_opt, ~, ~] = ...
    lambertMR(rE_dep.', rM_arr.', tof_opt_sec, muSun, orbitType, Nrev, Ncase, optionsLMR);

if ERRORopt ~= 0
    error('Lambert failed for optimal grid point (unexpected).');
end

%% Orbits of Earth and Mars (1 revolution each) using Kepler sampling
N_orb_pts = 400;
f_grid = linspace(0, 2*pi, N_orb_pts);

R_Earth = zeros(N_orb_pts,3);
R_Mars  = zeros(N_orb_pts,3);

for k = 1:N_orb_pts
    kepE_tmp       = kepE_dep;
    kepE_tmp(6)    = f_grid(k);        % true anomaly
    [rE_k, ~]      = kep2car_local(kepE_tmp, muSun);
    R_Earth(k,:)   = rE_k.';

    kepM_tmp       = kepM_arr;
    kepM_tmp(6)    = f_grid(k);
    [rM_k, ~]      = kep2car_local(kepM_tmp, muSun);
    R_Mars(k,:)    = rM_k.';
end

%% Transfer arc propagation (two-body about the Sun)
tspanT = linspace(0, tof_opt_sec, 400);
y0_T   = [rE_dep; v1_t_opt.'];

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, YT] = ode113(@(t,y) twoBodyODE_local(t,y,muSun), tspanT, y0_T, opts);
R_trans = YT(:,1:3);

%% Plot trajectory
figure; hold on; grid on; axis equal;
plot3(R_Earth(:,1), R_Earth(:,2), R_Earth(:,3), 'b', 'LineWidth', 1.5);  % Earth orbit
plot3(R_Mars(:,1),  R_Mars(:,2),  R_Mars(:,3),  'r', 'LineWidth', 1.5);  % Mars orbit
plot3(R_trans(:,1), R_trans(:,2), R_trans(:,3), 'g', 'LineWidth', 2.0);  % transfer arc

% Mark positions
plot3(rE_dep(1), rE_dep(2), rE_dep(3), 'bo', 'MarkerFaceColor','b', 'MarkerSize',7); % Earth @ dep
plot3(rM_arr(1), rM_arr(2), rM_arr(3), 'rs', 'MarkerFaceColor','r', 'MarkerSize',7); % Mars @ arr

% Sun at origin
plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',8);

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title(sprintf('Mars Express transfer (grid minimum) – \\DeltaV = %.3f km/s', dv_min));
legend('Earth orbit','Mars orbit','Transfer arc','Earth dep','Mars arr','Sun',...
       'Location','bestoutside');
view(3);
hold off;

%% ===== LOCAL FUNCTIONS =====

function mjd2000 = date2mjd2000_local(y, m, d, hh, mm, ss)
    % Convert calendar date/time to MJD2000 (days since 2000-01-01 12:00)
    if m <= 2
        y = y - 1;
        m = m + 12;
    end
    A  = floor(y/100);
    B  = 2 - A + floor(A/4);
    JD = floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + d + B - 1524.5;
    frac = (hh + mm/60 + ss/3600)/24;
    JD = JD + frac;
    mjd2000 = JD - 2451545.0;  % JD of 2000-01-01 12:00
end

function [r_I, v_I] = kep2car_local(kep, mu)
    % kep = [a e i RAAN omega theta] with angles in rad, a in km
    a    = kep(1);
    e    = kep(2);
    inc  = kep(3);
    RAAN = kep(4);
    om   = kep(5);
    f    = kep(6);

    p = a*(1 - e^2);

    % Perifocal coordinates
    r_pf = [ p*cos(f)/(1 + e*cos(f));
             p*sin(f)/(1 + e*cos(f));
             0 ];
    v_pf = [ -sqrt(mu/p)*sin(f);
              sqrt(mu/p)*(e + cos(f));
              0 ];

    % Rotation: perifocal -> inertial (ecliptic of date)
    R3_W = [ cos(RAAN) -sin(RAAN) 0;
             sin(RAAN)  cos(RAAN) 0;
             0          0         1];
    R1_i = [ 1     0          0;
             0 cos(inc) -sin(inc);
             0 sin(inc)  cos(inc)];
    R3_w = [ cos(om) -sin(om) 0;
             sin(om)  cos(om) 0;
             0        0       1];

    Q = R3_W * R1_i * R3_w;

    r_I = Q * r_pf;
    v_I = Q * v_pf;
end

function dy = twoBodyODE_local(~, y, mu)
    r = y(1:3);
    v = y(4:6);
    rnorm = norm(r);
    a = -mu * r / rnorm^3;
    dy = [v; a];
end