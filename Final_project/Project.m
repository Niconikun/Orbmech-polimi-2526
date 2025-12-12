clc
clear
close all

mu_Sun = astroConstants(4);

Earth = 3;
Mars = 4;
Jupiter = 5;
asteroid = 257323;

res = 10000;

Departure_time_earliest = [2030, 01, 01, 0, 0, 0.00];
Arrival_time_latest     = [2060, 01, 01, 0, 0, 0.00];

Departure_time_earliest = date2mjd2000(Departure_time_earliest);
Arrival_time_latest     = date2mjd2000(Arrival_time_latest);

Time_interval = linspace(Departure_time_earliest,Arrival_time_latest,res);

figure('Name','Solar system')
hold on
grid on
axis equal
ax = gca;
ax.Color = 'k';
ax.GridColor = [1 1 1];
view(3)

R_S = 10000000;
[X, Y, Z] = sphere(50);
X_S = R_S*X;
Y_S = R_S*Y;
Z_S = R_S*Z;

surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');

for i = 1:res

    time = Time_interval(i);
    
    title(sprintf('%.2f', mjd20002date(time)));

    kep_Earth    = uplanet(time, Earth);
    kep_Mars     = uplanet(time, Mars);
    kep_Jupiter  = uplanet(time, Jupiter);
    kep_asteroid = ephAsteroids(time, asteroid);

    plot_for_animated(kep_Earth,mu_Sun,1,kep_Earth(6),kep_Earth(6),'.b' )
    plot_for_animated(kep_Mars,mu_Sun,1,kep_Mars(6),kep_Mars(6),'.r' )
    plot_for_animated(kep_Jupiter,mu_Sun,1,kep_Jupiter(6),kep_Jupiter(6),'.m' )
    plot_for_animated(kep_asteroid,mu_Sun,1,kep_asteroid(6),kep_asteroid(6),'.' )

    pause(0.1);


end

%%

figure()
plot(1,1)


%% FUNCTIONS

function plot_for_animated(kep, mu, npoints, theta_i, theta_f, string)

theta = linspace(theta_i,theta_f,npoints);

R = zeros(3,npoints);

for i = 1:npoints

    kep(6) = theta(i);

    r = kep2car(kep,mu);

    R(:,i) = r;

end

plot3(R(1,:),R(2,:),R(3,:),string,'LineWidth',1)

end




 