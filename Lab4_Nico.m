clc
clear all
close all

%% Constants

deg2rad = @(x) x*pi/180;
day2sec = @(d) d*86400;


addpath("ephemerides\")
addpath("Functions\")
addpath("Functions\timeConversion\")

AU = 149597870.691;              % Astronomical Unit [km]


% Sun
mu_Sun = astroConstants(4);      % Sun's gravitational parameter [km^3/s^2];
R_Sun = astroConstants(3);       % Sun's radius [km]
M_Sun = mu_Sun / astroConstants(1);

% Earth
mu_Earth = astroConstants(13);      % Sun's gravitational parameter [km^3/s^2];
R_Earth = astroConstants(23);       % Sun's radius [km]
M_Earth = mu_Earth / astroConstants(1);
v_Earth = 29.78; %km/s            % Earth orbital velocity around the Sun [km/s]
v_inf_Earth = 11.186; %km/s        % Earth escape velocity [km/s]
r_SOI_Earth = SOI(M_Sun,M_Earth,AU);
theta_SOI_Earth = acos( (a*(1 - e^2)-r_SOI_Earth)/(r_SOI_Earth*e) );



% Mars
mu_Mars = astroConstants(14);      % Sun's gravitational parameter [km^3/s^2];
R_Mars = astroConstants(24);       % Sun's radius [km]
M_Mars = mu_Mars / astroConstants(1);
v_Mars = 24.07; %km/s              % Mars orbital velocity around the Sun [km/s]
v_inf_Mars = 5.027; %km/s          % Mars escape velocity [km/s]
r_SOI_Mars = SOI(M_Sun,M_Earth,AU);
theta_SOI_Mars = acos( (a*(1 - e^2)-r_SOI_Earth)/(r_SOI_Earth*e) );


%% Part 1a
v_inf_neg = [15.1;0;0];  %km/s
Delta = 9200;  %km
r_E = [1;0;0]*AU;  %km


% 1

[a, delta, e, rp] = solve_hyperbola_2d(norm(v_inf_neg),Delta,mu_Earth);


% 2
% In front of the planet
u = [0; 0; -1];
v_front = rodrigues_rotation(u,v_inf_neg,delta);

% Behind of the planet
u = [0; 0; 1];
v_behind = rodrigues_rotation(u,v_inf_neg,delta);
% Under of the planet
u = [0; -1; 0];
v_under = rodrigues_rotation(u,v_inf_neg,delta);

% 3

V_T = sqrt(mu_Sun / AU) * [0; 1; 0];

V_IN = V_T + v_inf_neg;

% In front of the planet
V_FRONT = V_T + v_front;

% Behind of the planet
V_BEHIND = V_T + v_behind;

% Under of the planet
V_UNDER = V_T + v_under;

% 4

G = astroConstants(2);
M_Sun = mu_Sun / G;
M_Earth = mu_Earth / G;
R_T = astroConstants(23);
R_S = 10000000;

r_SOI_Earth = SOI(M_Sun,M_Earth,AU);
theta_SOI = acos( (a*(1 - e^2)-r_SOI_Earth)/(r_SOI_Earth*e) );


% In front of the planet
[i_front, Omega_front, omega_front] = plane(v_inf,v_front,delta);
kep_front = [a,e,i_front,Omega_front,omega_front, 0];
KEP_FRONT_BEFORE = car2kep_mod(r,V_IN,mu_Sun);
KEP_FRONT_AFTER  = car2kep_mod(r,V_FRONT,mu_Sun);


% Behind the planet
[i_behind, Omega_behind, omega_behind] = plane(v_inf,v_behind,delta);
kep_behind = [a,e,i_behind,Omega_behind,omega_behind, 0];
KEP_BEHIND_BEFORE = car2kep_mod(r,V_IN,mu_Sun);
KEP_BEHIND_AFTER  = car2kep_mod(r,V_BEHIND,mu_Sun);

% Under of the planet
[i_under, Omega_under, omega_under] = plane(v_inf,v_under,delta);
kep_under = [a,e,i_under,Omega_under,omega_under, 0];
KEP_UNDER_BEFORE = car2kep_mod(r,V_IN,mu_Sun);
KEP_UNDER_AFTER  = car2kep_mod(r,V_UNDER,mu_Sun);


[X, Y, Z] = sphere(50);
X_T = R_T*X;
Y_T = R_T*Y;
Z_T = R_T*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

X_S = R_S*X;
Y_S = R_S*Y;
Z_S = R_S*Z;


% Plot in front of the planet 
figure('Name','In front of the planet')

subplot(1,2,1)
title('Heliocentric frame')
hold on
surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');
my_plotOrbit(KEP_FRONT_BEFORE, mu_Sun, 1000, KEP_FRONT_BEFORE(6)-3/2*pi, KEP_FRONT_BEFORE(6)      , 'b')
my_plotOrbit(KEP_FRONT_AFTER , mu_Sun, 1000, KEP_FRONT_AFTER(6)        , 2*pi - KEP_FRONT_AFTER(6), 'r')
my_plotOrbit(KEP_FRONT_BEFORE, mu_Sun, 1   , KEP_FRONT_BEFORE(6)       , KEP_FRONT_BEFORE(6)      , 'ok')
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

subplot(1,2,2)
title('Planetocentric frame')
hold on
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
my_plotOrbit(kep_front, mu_Earth, 1000, -1.5, +1.5, 'b')
my_plotOrbit(kep_front, mu_Earth, 1, 0, 0, 'ob')
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Plot behind the planet 
figure('Name','Behind the planet')

subplot(1,2,1)
title('Heliocentric frame')
hold on
surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');
my_plotOrbit(KEP_BEHIND_BEFORE, mu_Sun, 1000, KEP_BEHIND_BEFORE(6)-3/2*pi, KEP_BEHIND_BEFORE(6)      , 'b')
my_plotOrbit(KEP_BEHIND_AFTER , mu_Sun, 1000, KEP_BEHIND_AFTER(6)        , pi-KEP_BEHIND_AFTER(6), 'r')
my_plotOrbit(KEP_BEHIND_BEFORE, mu_Sun, 1   , KEP_BEHIND_BEFORE(6)       , KEP_BEHIND_BEFORE(6)      , 'ok')
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

subplot(1,2,2)
title('Planetocentric frame')
hold on
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
my_plotOrbit(kep_behind, mu_Earth, 1000, -1.5, +1.5, 'b')
my_plotOrbit(kep_behind, mu_Earth, 1, 0, 0, 'ob')
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Plot under the planet 
figure('Name','Under the planet')

subplot(1,2,1)
title('Heliocentric frame')
hold on
surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');
my_plotOrbit(KEP_UNDER_BEFORE, mu_Sun, 1000, KEP_UNDER_BEFORE(6)-3/2*pi, KEP_UNDER_BEFORE(6)      , 'b')
my_plotOrbit(KEP_UNDER_AFTER , mu_Sun, 1000, KEP_UNDER_AFTER(6)        , pi                       , 'r')
my_plotOrbit(KEP_UNDER_BEFORE, mu_Sun, 1   , KEP_UNDER_BEFORE(6)       , KEP_UNDER_BEFORE(6)      , 'ok')
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(3)

subplot(1,2,2)
title('Planetocentric frame')
hold on
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
my_plotOrbit(kep_under, mu_Earth, 1000, -1.5, +1.5, 'b')
my_plotOrbit(kep_under, mu_Earth, 1, 0, 0, 'ob')
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(3)


%% Function: solve_hyperbola_2d
% Solves for the parameters of a hyperbolic trajectory in 2D for flyby maneuvers.
% This function calculates the semi-major axis, turn angle, eccentricity, and periapsis radius
% of a hyperbolic orbit given the incoming velocity at infinity, the turn angle, and the gravitational parameter.
%
% Inputs:
%   v_inf_meno - Incoming velocity at infinity [km/s]
%   Delta      - Turn angle of the hyperbolic trajectory [rad]
%   mu         - Gravitational parameter of the central body [km^3/s^2]
%
% Outputs:
%   a          - Semi-major axis of the hyperbola [km] (negative for hyperbolas)
%   delta      - Turn angle [rad]
%   e          - Eccentricity of the hyperbola
%   rp         - Periapsis radius [km]
%
% Related to interplanetary flybys: Used to compute the geometry of the hyperbolic approach and departure trajectories
% during a planetary flyby, where the spacecraft's path is deflected by the planet's gravity.
function [a, delta, e, rp] = solve_hyperbola_2d(v_inf_meno,Delta,mu)

    a = -mu / v_inf_meno^2;

    delta = 2 * atan(-a / Delta );

    e = 1 / sin(delta / 2);

    rp = - mu / v_inf_meno^2 * ( 1 - e );

end




%% Function: rodrigues_rotation
% Performs a rotation of a vector using Rodrigues' rotation formula.
% This is used to rotate vectors in 3D space around a given axis by a specified angle.
%
% Inputs:
%   u     - Unit vector representing the axis of rotation
%   v     - Vector to be rotated
%   delta - Angle of rotation [rad]
%
% Outputs:
%   v_rot - Rotated vector
%
% Related to interplanetary exploration: Useful for transforming velocity vectors during flyby calculations,
% such as rotating the incoming and outgoing velocity vectors in the hyperbolic trajectory.
function v_rot = rodrigues_rotation(u,v,delta)

    v_rot = v*cos(delta) + cross(u,v)*sin(delta) + u*dot(u,v)*(1-cos(delta));

end



%% Function: my_plotOrbit
% Plots a 3D orbit trajectory given Keplerian elements over a range of true anomalies.
% This function generates points along the orbit and plots them in 3D space.
%
% Inputs:
%   kep     - Keplerian elements vector [a, e, i, Omega, omega, theta] where theta is initial true anomaly [rad]
%   mu      - Gravitational parameter [km^3/s^2]
%   npoints - Number of points to plot
%   theta_i - Initial true anomaly [rad]
%   theta_f - Final true anomaly [rad]
%   string  - Plot style string (e.g., 'b-' for blue line)
%
% Outputs:
%   None (plots directly)
%
% Related to interplanetary exploration: Used to visualize spacecraft orbits, including hyperbolic trajectories
% during flybys or elliptical orbits in interplanetary transfers.
function my_plotOrbit(kep, mu, npoints, theta_i, theta_f, string)

    theta = linspace(theta_i,theta_f,npoints);

    R = zeros(3,npoints);

    for i = 1:npoints

        kep(6) = theta(i);

        r = kep2car(kep,mu);

        R(:,i) = r;

    end

    plot3(R(1,:),R(2,:),R(3,:),string,'LineWidth',1.5)

end


%% Function: SOI
% Calculates the radius of the Sphere of Influence (SOI) of a smaller body in the gravitational field of a larger one.
% The SOI defines the region where the smaller body's gravity dominates over the larger one.
%
% Inputs:
%   M        - Mass of the larger body (e.g., Sun) [kg]
%   m        - Mass of the smaller body (e.g., planet) [kg]
%   distance - Distance between the two bodies [km]
%
% Outputs:
%   r_SOI    - Radius of the Sphere of Influence [km]
%
% Related to interplanetary exploration: Determines the boundary beyond which the planet's gravity
% can be neglected for trajectory calculations, crucial for patched conics approximations in flybys.
function [r_SOI] = SOI(M,m,distance)
    r_SOI = distance*(m/M)^(2/5);
end




%% Function: plane
% Determines the orbital plane parameters (inclination, right ascension of ascending node, argument of periapsis)
% from the incoming and outgoing velocity vectors at infinity and the turn angle.
% This is used to define the plane of the hyperbolic trajectory during a flyby.
%
% Inputs:
%   v_inf_meno - Incoming velocity vector at infinity [km/s]
%   v_inf_piu  - Outgoing velocity vector at infinity [km/s]
%   delta      - Turn angle of the hyperbolic trajectory [rad]
%
% Outputs:
%   i      - Inclination of the orbital plane [rad]
%   Omega  - Right ascension of the ascending node [rad]
%   omega  - Argument of periapsis [rad]
%
% Related to interplanetary exploration: Essential for calculating the orientation of the flyby trajectory
% relative to the ecliptic or other reference planes, aiding in mission planning for planetary encounters.
function [i, Omega, omega] = plane(v_inf_meno,v_inf_piu,delta)

    toll = 1e-24;

    x_vec = [ 1; 0; 0 ];
    z_vec = [ 0; 0; 1 ];

    h_dir = cross(v_inf_meno,v_inf_piu) / norm(cross(v_inf_meno,v_inf_piu));

    i = acos(h_dir(3));

    N = cross(z_vec, h_dir);
    N_norm = norm(N);
    if N_norm < toll
        N = x_vec;
        N_norm = 1;
    end

    if isequal(N, x_vec)
        Omega = 0;
    else
        if N(2) >= 0
            Omega = acos(N(1) / N_norm);
        else
            Omega = 2*pi - acos(N(1) / N_norm);
        end
    end

    beta = ( pi - delta ) / 2;

    v_inf_piu_norm = v_inf_piu / norm(v_inf_piu);

    e_vec_norm = - rodrigues_rotation(h_dir,v_inf_piu_norm,beta);

    if (i ~= 0) && (i ~= pi)
        if e_vec_norm(3) >= 0
            omega = acos( dot(N, e_vec_norm) / (N_norm) );
        else
            omega = 2*pi - acos( dot(N, e_vec_norm) / (N_norm) );
        end

    elseif i == 0
        if e_vec_norm(2) >= 0
            omega = acos( dot(N, e_vec_norm) / (N_norm) );
        else
            omega = 2*pi - acos( dot(N, e_vec_norm) / (N_norm) );
        end

    elseif i == pi
        if e_vec_norm(2) <= 0
            omega = acos( dot(N, e_vec_norm) / (N_norm) );
        else
            omega = 2*pi - acos( dot(N, e_vec_norm) / (N_norm) );
        end
    end



end