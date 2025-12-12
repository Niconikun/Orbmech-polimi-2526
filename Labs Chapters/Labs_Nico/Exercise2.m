%% Exercise 5

%the time provided is the winter solstice, thus the ascension is 1.5*pi and
%the declination is 23.4 deg.

delta = 23.4; %deg, declination
ascension = * pi; % radians, ascension
elevation = 0; %deg, elevation
h = 0;
phi = 45;
%calculate the sun's azimuth
azimuth = atan2d(sind(ascension), cosd(ascension) * sind(delta) - tand(elevation) * cosd(ascension));
% Display the calculated sun azimuth
fprintf('The sun azimuth is %.2f degrees.\n', azimuth);


%% Exercise 6


%% Exercise 8
telescope_latitude = 35.14; %deg
star_declination = 13.53; %deg
%calculate the mazimum elevation of the star above the horizon
max_elevation = asin(sin(deg2rad(star_declination)) * sin(deg2rad(telescope_latitude))) * (180/pi);
% Display the maximum elevation of the star above the horizon
fprintf('The maximum elevation of the star is %.2f degrees.\n', max_elevation);

%% Exercise 9

alpha = 189.25; %deg
delta = 9.5; %deg
epsilon = 23.44; %deg
beta = 12.3857; %deg

lambda = 90 - acosd((cosd(90-delta)-cosd(epsilon)*cosd(90-beta))/(sind(epsilon)*sind(90-beta)));
fprintf(lambda + "°");

%% Exercise 14

%% Exercise 15

phi = 55.75; %deg
lambda = 2.51; %hrs
delta_sun = 7.81; %deg
E_28_2_mins = -12;
E_28_2_sec = 27;
E_1_3_mins = -12;
E_1_3_sec = 16;

h_1 = 0;

Az = acosd(cosd(90-delta_sun)/sind(90-delta_sun));
t = 360- acosd(-cosd(90-delta_sun)*cosd(90-phi)/(sind(90-delta_sun)*sind(90-phi)));
fprintf('Azimuth for h_1: %.2f°, Hour angle for h_1: %.2f°\n', Az, t)

h_2 = 180; % Set the second hour angle
Az_2 = 360 - acosd(cosd(90 - delta_sun) / sind(90 - delta_sun)); % Calculate the azimuth for h_2
t_2 = acosd(-cosd(90 - delta_sun) * cosd(90 - phi) / (sind(90 - delta_sun) * sind(90 - phi))); % Calculate the hour angle for h_2
fprintf('Azimuth for h_2: %.2f°, Hour angle for h_2: %.2f°\n', Az_2, t_2);

%LAT and LCT


%% 2.16:
lat = 34*pi/180;
long = 160*pi/180;
long_h = 160/15;
LCT = 5 + 16/60 + 10/3600;
dec = (8 + 19/60)*pi/180;
E12 = - (13/3600);
E13 = 6/3600;
E14 = 25/3600;

GMT = LCT + ceil((long*180/pi - 7.5)/15); % at 12-04
LAT = GMT*((E13 - E14)/24 + 1) + long_h + E12;
ts_h = LAT - 12;
ts_deg = ts_h*15
ts = ts_deg*pi/180;

h = pi/2 - acos(cos(pi/2 + dec)*cos(pi/2 - lat) + sin(pi/2 + dec)*sin(pi/2 - lat)*cos(ts));
h_deg = h*180/pi

Az = pi + acos((cos(pi/2 + dec) - cos(pi/2 - h)*cos(pi/2 - lat))/(sin(pi/2 - h)*sin(pi/2 - lat)));
Az_deg = Az*180/pi


%% Exercise 17

