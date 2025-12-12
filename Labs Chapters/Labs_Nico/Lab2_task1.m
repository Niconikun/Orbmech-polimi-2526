clc
clearvars
addpath("Functions 2bp\");

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);
J2 = astroConstants(33);
% Initial condition

kep_parameters = [8350,0.1976,60,270,45,230]; % a [km],e [-],i [deg],Omega [deg],omega [deg],theta [deg]

y0 = kep2car(kep_parameters, mu_E)';

% Set time span
a = kep_parameters(1); % Semi-major axis [km]
T = orbital_period(mu_E,a); % Orbital period [s]
tspan = linspace( 0, 2*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bpp(t,y,mu_E, J2, R_E), tspan, y0, options );

r = Y(:,1:3);
v = Y(:,4:6);
textureImage = imread('EarthTexture.jpg');

% Create the sphere
[XS,YS,ZS] = sphere(50);

scaleFactor = 6371; % Radius of Earth in kilometers; adjust as needed

% Scale the sphere
XS = XS * scaleFactor;
YS = YS * scaleFactor;
ZS = ZS * -scaleFactor;

% Plot the results
figure(1)
axis equal;
clf;

% Earth
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none');
hold on;

% Ideal orbit with color  
for i = 1:length(T)-1
    plot3(Y(i:i+1,1), Y(i:i+1,2), Y(i:i+1,3), 'r');
end

hold off;

figure(2)
% get specificenergy from all variables of r and v
specificEnergy = zeros(size(tspan)); % Preallocate for specific energy values
for i = 1:length(tspan)
    specificEnergy(i) =+ specificenergy(norm(Y(i,4:6)), mu_E, norm(Y(i,1:3)), a);
end
plot(T,specificEnergy);
xlabel('X [km]'); ylabel('Y [km]');
title('Specific Energy');
grid on;

figure(3)
[h_v,h,h_u] = angular_momentum(r,v);
hold on;
plot(T,h_v);
H = ones(length(tspan)).*h;
plot(T,H);
hold off;

xlabel('X [km]'); ylabel('Y [km]');
title('Angular Momentum');
grid on;

figure(4)
[e_v, e, e_unit] = eccentricity(r,v,h_v,mu_E);
plot(T,norm(e_v(:,1)),T,norm(e_v(:,2)),T,norm(e_v(:,3)));
hold on;
e_module = ones(length(tspan)).*e;
plot(T,e_module)
xlabel('X [km]'); ylabel('Y [km]');
title('Eccentricity Vector');
grid on;

figure(5)
%get the velocity components using getvelcomponents for each vector in time
v_radial = ones(length(tspan));
v_transversal = ones(length(tspan),3);% Preallocate for velocity components
for i = 1:length(tspan)
    [radial_vector,unit_radial_vector, transversal_vector, unit_transversal_vector] = getvelcomponents(h_v(i,:), v(i,:), r(i,:));
    v_radial(i,:) =+ norm(radial_vector);
    v_transversal(i,:) = norm(transversal_vector);
end

plot(T,v_transversal,T,v_radial);

xlabel('X [km]'); ylabel('Y [km]');
title('Velocity components');
grid on;

figure(6)
dotproduct = zeros(length(tspan));
for i=1:length(tspan)
    dotproduct(i) =+ dot(h_v(i,:),e_v(i,:));
end
plot(T, dotproduct);
xlabel('X [km]'); ylabel('Y [km]');
title('e-h dot product');
grid on;
