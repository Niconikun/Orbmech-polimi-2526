function [alpha, delta, lon, lat] = groundTrack(r, time, w_planet, t_0, theta_g_t_0)
%GROUNDTRACK Computes the ground track of a satellite given its position vector.

    R = norm(r);
    fprintf('Position magnitude: %f km\n', R); % Debug
    
    % Calculate declination (latitude)
    delta = asin(r(3) / R);
    fprintf('Declination (delta): %f rad, %f deg\n', delta, rad2deg(delta));
    
    % Calculate right ascension
    alpha = atan2(r(2), r(1));
    fprintf('Right ascension (alpha): %f rad, %f deg\n', alpha, rad2deg(alpha));
    
    % Calculate Greenwich sidereal time
    w_planet_rad = deg2rad(w_planet);
    theta_g = theta_g_t_0 + w_planet_rad * (time - t_0);
    fprintf('Greenwich sidereal time (theta_g): %f rad, %f deg\n', theta_g, rad2deg(theta_g));
    
    % Calculate longitude
    lon = rad2deg(alpha - theta_g);
    lon = mod(lon + 180, 360) - 180;
    
    % Calculate latitude
    lat = rad2deg(delta);
    
    fprintf('Final output - Lon: %f deg, Lat: %f deg\n\n', lon, lat);
end