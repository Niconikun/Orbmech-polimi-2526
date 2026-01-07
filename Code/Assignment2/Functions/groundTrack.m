function [alpha, delta, lon, lat] = groundTrack(S, lon_G0, omega_E, t_vec, str)
% groundTrack Ground Track of a Keplerian orbit
% 
% PROTOTYPE
%   [alpha, delta, lon, lat] = groundTrack(S, lon_G0, omega_E, t_vec, str)
%
% INPUT
%   S [Nx3]          Vector of position
%   lon_G0 [1]       Longitude of the Greenwich meridian at initial time
%   omega_E [1]      Angular velocity of the planet
%   t_vec [1xN]      Vector of instant of time
%   str [1xM]        String for Color and Marker definition for the plot
% 
% OUTPUT
%   alpha [-]        Right ascensions in Earth Centred Equatorial Inertial frame
%   delta [-]        Declination in Earth Centred Equatorial Inertial frame
%   lon [-]          Longitude with respect to rotating Earth
%   lat [-]          Latitude with respect to rotating Earth
%
% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-10-21: First version
%   2025-10-24: Second version
%   2025-11-13: Third version
%   2025-11-15: Fourth version
% 
% -------------------------------------------------------------------------

% 1. Check
if isempty(t_vec)
        error('groundTrack: t_vec is empty');
end
if ~isvector(t_vec)
    error('groundTrack: t_vec is not a vector');
end
if size(S,1) ~= numel(t_vec)
    error('groundTrack: number of row of S (%d) is different from the length of t_vec (%d)', size(S,1), numel(t_vec));
end

[~, cols] = size(S);
if cols < 3
    error('S must have at least 3 columns (x, y, z)');
end

% 2. Position
r = S(:, 1:3);

% 3. Normalization 
r_norm = sqrt(sum(r.^2, 2));

% 4. Declination
delta = asin(r(:, 3) ./ r_norm)';

% 5. Latitude 
lat = delta;

% 6. Right ascension 
alpha = atan2(r(:, 2), r(:, 1))';
alpha = wrapTo2Pi(alpha);

% 7. Longitude
lon_G = zeros(1, size(r, 1));
for i = 1:size(r, 1)
    lon_G(i) = (omega_E)*t_vec(i) + lon_G0;
end
lon = alpha - lon_G; 
lon = wrapToPi(lon);
lon_periodic = wrapToPi(lon);

% 8. Adjusting longitude for periodicity
%lon_periodic = mod(lon + pi, 2*pi) - pi; % Wraps longitude to [-pi, pi]

% 9. Plot of the ground Track
groundTrack_plot(lon_periodic, lat, str)

end

function groundTrack_plot(lon, lat, str)

    % 1. Definition of Marker and Color
    valid_colors = 'rgbcmyk';  
    valid_markers = '+o*.xsd^v><ph'; 

    % 2. Take marker and color
    if nargin >= 3 && ~isempty(str)

        c = str(end);
        if contains(valid_colors, c)
            trackColor = c;
        else
            trackColor = 'r';
        end

        m = str(1);
        if contains(valid_markers, m)
            trackMarker = m;
        else
            trackMarker = '.';
        end

    else

        trackColor = 'r';
        trackMarker = '.';

    end

    % 3, Transformation of latitude and longitudine in [deg]
    x_plot = lon * 180/pi;
    y_plot = lat * 180/pi;

    % 4. Plot
    hold on
    plot(x_plot, y_plot, trackMarker, 'LineStyle', 'none', 'Color', trackColor)
    plot(x_plot(1), y_plot(1), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', trackColor)
    text(x_plot(1) + 3, y_plot(1), 'START', 'FontWeight', 'bold', 'Color', trackColor)
    plot(x_plot(end), y_plot(end), 's', 'MarkerSize', 10, 'LineWidth', 2, 'Color', trackColor)
    text(x_plot(end) + 3, y_plot(end), 'END', 'FontWeight', 'bold', 'Color', trackColor)
    grid on
    axis on
    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    axis([-180 180 -90 90]);
    set(gca, 'YDir', 'normal');

end