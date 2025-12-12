%% Exercise 3.1 
clc
clear all

mu_E = 398600; %km3/s2
T_GEO = 23*60*60 + 56*60 + 4;

omega_E = 2*pi/T_GEO; %rad/s;

r_GEO = 42164; %km
v_GEO = 3.0747; %km/s

% the period of the phasing orbit can be determined like this

T2 = (20*pi/180 + 6*pi)/(3*omega_E); %s

a2 = (T2*sqrt(mu_E)/(2*pi)).^(2/3);

r_c = 2 * a2 - r_GEO;

h2 = sqrt(2*mu_E)*sqrt(r_GEO*r_c/(r_GEO + r_c));

v2 = h2/r_GEO;

deltav = 2*(v2-v_GEO);

fprintf('Delta V: %f', deltav);

%% Exercise 3.4
clc
clear all
e1 = 0.1755;
e2 = e1;
r_pericenter = 1400 + 6378.316; %km
mu_E = 398600; %km3/s2
alpha = 46.13;

theta_1 = alpha;
theta_2 = 0;



deltha_w = theta_2 - theta_1;

theta = deltha_w/2; % true anomaly

theta1 = 360 + theta; % true anomaly

theta2 = theta1 - 180; % true anomaly

fprintf("point 1: %f \n", theta1);
fprintf("point 2: %f \n", theta2);

    % Conversion to radians

theta1_rad = deg2rad(theta1); 
theta2_rad = deg2rad(theta2);

a = r_pericenter/(1-e1);

p = a*(1-e1^2);
 
deltha_V = 2*sqrt(mu_E/p)*e1*sin(theta2_rad);

E = atan(  tan(theta2_rad/2) * sqrt((1-e1) / (1+e1)) ) * 2; % Eccentric anomaly

M = E - e1*sin(E); % Mean anomaly

n = sqrt(mu_E/a^3); % Mean motion

Deltha_t = M /n;

Deltha_t_min = Deltha_t/60;

fprintf("time: %f [min] \n",Deltha_t_min);



%% Exercise 3.6
clc
clear all
mu_E = 398600; %km3/s2
angle_manoeuvre = 53; %°
deltav = 1.8; %km/s

r_apoapsis = 17000; %km
r_periapsis = 7000; %km

e1 = (r_apoapsis-r_periapsis)/(r_apoapsis+r_periapsis);

h1 = sqrt(r_periapsis*((1+e1*cosd(0))*mu_E));

v_vertical_1 = h1/r_periapsis;
v_radial_1 = 0;

deltav_vertical = deltav*cosd(angle_manoeuvre);
deltav_radial = deltav*sind(angle_manoeuvre);

h2 = h1 + r_periapsis*deltav_vertical;
A = (v_vertical_1+deltav_vertical)*(v_radial_1+deltav_radial)*v_vertical_1^2;
B = ((mu_E/r_periapsis)*(deltav_vertical*(2*v_vertical_1+deltav_vertical)+e1*cosd(0)*(v_vertical_1+deltav_vertical)^2));
C = A/B;
theta_2 = -atand(C);

%print theta_2
fprintf("point 2: %f \n", theta_2);
e2 = ((h1 + r_periapsis*deltav_vertical)^2*e1*cosd(0) + (2*h1 + r_periapsis*deltav_vertical)*r_periapsis*deltav_vertical)/(h1^2*cosd(theta_2));

r_periapsis_2 = h2^2 / (mu_E*(1+e2));
r_apoapsis_2 = h2^2 / (mu_E*(1-e2));

%print periapsis and apoapsis
fprintf('Periapsis 2: %f km\n', r_periapsis_2);
fprintf('Apoapsis 2: %f km\n', r_apoapsis_2);




%% Exercise 3.8
clc
clear all

R_E =6378.136; %km
mu_E = 398600; %km3/s2

r_p1 = 600 + R_E;
r_a1 = 800 + R_E;

e1 = (r_a1-r_p1)/(r_a1+r_p1);

h1 = sqrt(r_p1*((1+e1*cosd(0))*mu_E));

v_vertical_1 =  h1/r_p1;

r_p2 = r_p1;

r_a2 = R_E + 15000;

e2 = (r_a2 - r_p2) / (r_a2 + r_p2);
h2 = sqrt(r_p2 * ((1 + e2 * cosd(0)) * mu_E));
v_vertical_2 = h2 / r_p2;

fprintf('deltav1: %f \n', abs(v_vertical_2-v_vertical_1));

e3 = 0;

h3 = sqrt(r_a2 * mu_E);

v_vertical_3 = h3 / r_a2;

fprintf('deltav2: %f \n', abs(v_vertical_3-v_vertical_2));

%% Exercise 3.10
clc 
clear all
R_p1 = 7500; %km
mu_E = 398600; %km3/s2
R_p3 = 120000; %km
R_a2 = 190000; %km


h1 = sqrt(R_p1 * mu_E);

v_vertical_1 = h1 / R_p1;
e2 = (R_a2 - R_p1) / (R_a2 + R_p1);
h2 = sqrt(R_a2 * ((1 + e2 * cosd(180)) * mu_E));

v_vertical_2 = h2 / R_p1;

a_2 = (R_a2 + R_p1)/2;

ToF1 = 0.5*(2*pi*a_2^(1.5)/sqrt(mu_E));

deltav_1 = v_vertical_2 - v_vertical_1;


e3 = (R_a2 - R_p3) / (R_a2 + R_p3);

h3 = sqrt(R_a2 * ((1 + e3 * cosd(180)) * mu_E));
v_vertical_3 = h3 / R_a2;
deltav_2 = v_vertical_3 - v_vertical_2;

a_3 = (R_a2 + R_p3)/2;

ToF2 = 0.5*(2*pi*a_3^(1.5)/sqrt(mu_E));

h4 = sqrt(R_p3 * mu_E);
v_vertical_4 = h4 / R_p3;
deltav_3 = v_vertical_4 - v_vertical_3;

deltav_total = abs(deltav_1)+abs(deltav_2)+abs(deltav_3);

ToF_total = ToF1 + ToF2;


fprintf('Total Delta V: %f\n', deltav_total);
fprintf('Total Time of Flight: %f\n [days]', ToF_total/(60*60*24));

% Compare with a hohmann transfer


e2_h = (R_p3 - R_p1) / (R_p3 + R_p1);
h2_h = sqrt(R_p3 * ((1 + e2_h * cosd(180)) * mu_E));

v_vertical_2_h = h2_h / R_p1;

a_2_h = (R_p3 + R_p1)/2;

ToF1_h = 0.5*(2*pi*a_2_h^(1.5)/sqrt(mu_E));

deltav_1_h = v_vertical_2_h - v_vertical_1;


e3 = (R_a2 - R_p3) / (R_a2 + R_p3);

h3_h = sqrt(R_p3 * mu_E);
v_vertical_3_h = h3_h / R_p3;
deltav_2_h = v_vertical_3_h - v_vertical_2_h;

ToF_total_h = ToF1_h;

deltav_total_h = abs(deltav_1_h)+abs(deltav_2_h);

fprintf('Total Delta V (Hohmann): %f\n', deltav_total_h);
fprintf('Total Time of Flight (Hohmann): %f\n', ToF_total_h);

fprintf('Total Delta V Ratio (Hohmann/Bielliptic): %f\n', deltav_total/deltav_total_h);
fprintf('Total Time of Flight Ratio (Hohmann/Bielliptic): %f\n', ToF_total/ToF_total_h);






