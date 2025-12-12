clc
clear
close all

lon = linspace(-pi,pi,1000);
lat = sin(lon)*pi/4;

img = imread('EarthTexture.jpg');
img = flipud(img);
[y_length, x_length, ~] = size(img);

x_plot = x_length/2 + lon / pi * x_length / 2;
y_plot = y_length/2 - lat / pi * y_length;

figure
imshow(img, 'XData', [-180 180], 'YData', [-90 90]);
set(gca, 'YDir', 'normal');
hold on
plot(x_plot,y_plot,'.r','LineStyle','none')
grid on
axis on



%%
clc
clear
close all

exampleimage = imread('peppers.png');
exampleimage = flipud(exampleimage);
imshow(exampleimage)
hold on
x=1:20:300;
y=x;
scatter(x,y,'filled','w')
set(gca,'YDir','normal')


%% Plotting - From Nando

lambda = rad2deg(lambda);
% Adjusting longitude for periodicity
lambda_periodic = mod(lambda + 180, 360) - 180; % Wraps longitude to [-180, 180]

% Creating a 2D matrix with adjusted longitude and corresponding latitude
ground_track = [lambda_periodic', phi'];

figure()
plot(lambda_periodic, rad2deg(phi)) % Plotting the ground track
xlabel('Longitude [deg]'); ylabel('Latitude [deg]')
title('Spacecraft Ground Track');
grid on;