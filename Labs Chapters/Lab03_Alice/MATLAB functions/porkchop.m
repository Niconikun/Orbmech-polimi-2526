function [DV_TOT, minDV, best_date] = porkchop(dep, arr, day_step, planets, levels, v_inf)
% porkchop is a function that find the minimum delta v cost with a grid
% search algorithm for a grid of departure and arrival times covering the 
% time windows provided
% 
% PROTOTYPE
%  [DV_TOT, best_dep, best_arr, minDV, deps, arrs, DV1] = porkchop(dep, arr, day_step, planets, levels)
%
% INPUT:
%   dep [2x3]        Windows for departure in calendar date
%   arr [2x3]        Windows for arrivals in calendar date
%   day_step [1]     Step in days
%   planets [2x1]    Indexes of the two planets
%   levels [1xN]     Number of levels for the contour plot
%   v_inf [1]        Maximum excess velocity from launcher. Negative without considering launcher constraint
%
% OUTPUT:
%   DV_TOT [IxJ]     Matrix with the delta v cost
%   minDV [1]        Minimum delta v
%   best_date [2x6]  Dates of departure and arrival for the mission with minimum delta v
% 
% CONTRIBUTORS:
%   Alice Muscas
%
% VERSIONS
% 2025-11-23: First Version
% 
% -------------------------------------------------------------------------

% Sun planetary constant [km^3/s^2]
muS = astroConstants(4);

% Extraction of Data
dep_1 = dep(1, :);
dep_2 = dep(2, :);
arr_1 = arr(1, :);
arr_2 = arr(2, :);
planet_1 = planets(1);
planet_2 = planets(2); 

% Conversion of date in MJD2000
mjd_dep_1 = date2mjd2000(dep_1);
mjd_dep_2 = date2mjd2000(dep_2);
mjd_arr_1 = date2mjd2000(arr_1);
mjd_arr_2 = date2mjd2000(arr_2);

% Grid search for day
deps = mjd_dep_1 : 1 : mjd_dep_2;
arrs = mjd_arr_1 : 1 : mjd_arr_2;
DV_TOT = zeros(length(deps), length(arrs));
DV1 = zeros(length(deps), length(arrs));

for i = 1:length(deps)
    for j = 1:length(arrs)

        if arrs(j) < deps(i)
            DV_TOT(i, j) = NaN;
            DV1(i, j) = NaN;
        else
            % Keplerian element vector
            [kep_dep, ~] = uplanet(deps(i), planet_1);
            [r_dep, v_dep] = kep2car(kep_dep, muS);
            [kep_arr, ~] = uplanet(arrs(j), planet_2);
            [r_arr, v_arr] = kep2car(kep_arr, muS);
            
            % Time of flight in seconds
            ToF = (arrs(j) - deps(i)) * 3600 * 24; % [s]
    
            % Velocity cost
            [dv_tot, dv1, dv2] = VelocityCost(r_dep, v_dep, r_arr, v_arr, ToF, muS); % [km/s]
            DV_TOT(i, j) = dv_tot;
            DV1(i, j) = dv1;
        end

    end
end

% Plot
figure('Name', 'Porckchop plot of the Mission'); 
hold on;
[C, h] = contour(deps, arrs, DV_TOT', levels);
clabel(C, h);
colorbar;
xlabel('Departure MJD2000');
ylabel('Arrival MJD2000');
title('Porkchop Plot');

% Constrained case
if v_inf >= 0
    valid = DV1 <= v_inf;
    DV_TOT(~valid) = NaN;
end

% Minimum delta v
[minDV, idx] = min(DV_TOT(:));
if isnan(minDV)
    error("No available solution")
end
[i, j] = ind2sub(size(DV_TOT), idx);

% Dates for best solution
best_dep = deps(i);
best_arr = arrs(j);

% Grid search for hour
deps = best_dep - 1 : 1/24 : best_dep + 1;
arrs = best_arr - 1 : 1/24 : best_arr + 1;
DV_TOT = zeros(length(deps), length(arrs));
DV1 = zeros(length(deps), length(arrs));

for i = 1:length(deps)
    for j = 1:length(arrs)

        if arrs(j) < deps(i)
            DV_TOT(i, j) = NaN;
            DV1(i, j) = NaN;
        else
            % Keplerian element vector
            [kep_dep, ~] = uplanet(deps(i), planet_1);
            [r_dep, v_dep] = kep2car(kep_dep, muS);
            [kep_arr, ~] = uplanet(arrs(j), planet_2);
            [r_arr, v_arr] = kep2car(kep_arr, muS);
            
            % Time of flight in seconds
            ToF = (arrs(j) - deps(i)) * 3600 * 24; % [s]
    
            % Velocity cost
            [dv_tot, dv1, dv2] = VelocityCost(r_dep, v_dep, r_arr, v_arr, ToF, muS); % [km/s]
            DV_TOT(i, j) = dv_tot;
            DV1(i, j) = dv1;
        end

    end
end

% Constrained case
if v_inf >= 0
    valid = DV1 <= v_inf;
    DV_TOT(~valid) = NaN;
end

% Minimum delta v
[minDV, idx] = min(DV_TOT(:));
if isnan(minDV)
    error("No available solution")
end
[i, j] = ind2sub(size(DV_TOT), idx);

% Dates for the new best solution
best_dep = deps(i);
best_arr = arrs(j);

% Grid search for minutes
deps = best_dep - 1/24 : 1/24/60/4 : best_dep + 1/24;
arrs = best_arr - 1/24 : 1/24/60/4 : best_arr + 1/24;
DV_TOT = zeros(length(deps), length(arrs));
DV1 = zeros(length(deps), length(arrs));

for i = 1:length(deps)
    for j = 1:length(arrs)

        if arrs(j) < deps(i)
            DV_TOT(i, j) = NaN;
            DV1(i, j) = NaN;
        else
            % Keplerian element vector
            [kep_dep, ~] = uplanet(deps(i), planet_1);
            [r_dep, v_dep] = kep2car(kep_dep, muS);
            [kep_arr, ~] = uplanet(arrs(j), planet_2);
            [r_arr, v_arr] = kep2car(kep_arr, muS);

            % Time of flight in seconds
            ToF = (arrs(j) - deps(i)) * 3600 * 24; % [s]

            % Velocity cost
            [dv_tot, dv1, dv2] = VelocityCost(r_dep, v_dep, r_arr, v_arr, ToF, muS); % [km/s]
            DV_TOT(i, j) = dv_tot;
            DV1(i, j) = dv1;
        end

    end
end

% Constrained case
if v_inf >= 0
    valid = DV1 <= v_inf;
    DV_TOT(~valid) = NaN;
end

% Minimum delta v
[minDV, idx] = min(DV_TOT(:));
if isnan(minDV)
    error("No available solution")
end
[i, j] = ind2sub(size(DV_TOT), idx);

% Dates for the new best solution
best_dep = deps(i);
best_arr = arrs(j);


% Grid search for seconds
deps = best_dep - 1/24/60 : 1/24/3600/4 : best_dep + 1/24/60;
arrs = best_arr - 1/24/60 : 1/24/3600/4 : best_arr + 1/24/60;
DV_TOT = zeros(length(deps), length(arrs));
DV1 = zeros(length(deps), length(arrs));

for i = 1:length(deps)
    for j = 1:length(arrs)

        if arrs(j) < deps(i)
            DV_TOT(i, j) = NaN;
            DV1(i, j) = NaN;
        else
            % Keplerian element vector
            [kep_dep, ~] = uplanet(deps(i), planet_1);
            [r_dep, v_dep] = kep2car(kep_dep, muS);
            [kep_arr, ~] = uplanet(arrs(j), planet_2);
            [r_arr, v_arr] = kep2car(kep_arr, muS);

            % Time of flight in seconds
            ToF = (arrs(j) - deps(i)) * 3600 * 24; % [s]

            % Velocity cost
            [dv_tot, dv1, dv2] = VelocityCost(r_dep, v_dep, r_arr, v_arr, ToF, muS); % [km/s]
            DV_TOT(i, j) = dv_tot;
            DV1(i, j) = dv1;
        end

    end
end

% Constrained case
if v_inf >= 0
    valid = DV1 <= v_inf;
    DV_TOT(~valid) = NaN;
end

% Minimum delta v
[minDV, idx] = min(DV_TOT(:));
if isnan(minDV)
    error("No available solution")
end
[i, j] = ind2sub(size(DV_TOT), idx);

% Dates for the new best solution
best_dep = deps(i);
best_arr = arrs(j);


best_date = [best_dep; best_arr];

end