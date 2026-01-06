function theta_G = Greenwich_longitude(date)
% Greenwich_longitude calculates the Greenwich Mean Longitude for a given date.
% 
% Inputs:
%   date - A datetime object representing the date and time for which the
%          Greenwich longitude is to be calculated.
%
% Outputs:
%   theta_G - The calculated Greenwich Mean Longitude in degrees.
%
% The function computes the Julian date from the input date, converts it to
% Julian centuries, and then calculates the mean longitude based on the
% provided date and time. The result is adjusted for the time of day.

Y = year(date);
M = month(date);
D = day(date);
h = hour(date);
m = minute(date);
s = second(date);

J0 = 367*Y - floor(7 * (Y + floor((M + 9) / 12)) / 4) + floor(275 * M / 9) + D + 1721013.5;

T0 = (J0 - 2451545) / 36525;
theta_G0 = 100.4606184 + 36000.77004 * T0 + 0.000387933 * T0^2 - 2.583 * T0^3 * 1e-8;

theta_G0 = mod(theta_G0, 360);

theta_G = theta_G0 + 360.98564724*(h + m/60 + s/3600)/24; % Adjust for time of day

end
