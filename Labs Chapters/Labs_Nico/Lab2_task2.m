clc
clearvars
addpath("Functions 2bp\");

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(29);
J2_E = astroConstants(9);
% Initial condition
r0 = [ 6495; -970; -3622 ]; % [km]
v0 = [ 4.752; 2.130; 7.950]; % [km/s] 
y0 = [ r0; v0];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E); % Semi-major axis [km]
T_period = orbital_period(mu_E,a); % Orbital period [s]
tspan = linspace( 0, 50*T_period, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bpp(t,y,mu_E,J2_E,R_E), tspan, y0, options );
[ T2, Y2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

r = Y(:,1:3);
r2 = Y2(:,1:3);
v = Y(:,4:6);
v2 = Y2(:,4:6);
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

hold on
plot3( Y(:,1), Y(:,2), Y(:,3))
plot3( Y2(:,1), Y2(:,2), Y2(:,3), 'r-' )
surf(XS, YS, ZS, 'FaceColor', 'texturemap', 'CData', textureImage, 'EdgeColor', 'none');

ax = gca;
C = T2;
surface([Y2(:,1) Y2(:,1)], [Y2(:,2) Y2(:,2)], [Y2(:,3) Y2(:,3)], ...
        [C C], 'FaceColor','none', 'EdgeColor','interp', 'LineWidth',1.5);

clim(ax, [min(C) max(C)])

createcolorbar(ax)


% Finalize the plot
hold off;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;


figure(2)
% get specificenergy from all variables of r and v
specificEnergy1 = zeros(length(T), 1);
specificEnergy2 = zeros(length(T2), 1);
for i = 1:length(T)
    specificEnergy1(i) = specificenergy(norm(v(i,:)), mu_E, norm(r(i,:)), []);
end
for i = 1:length(T2)
    specificEnergy2(i) = specificenergy(norm(v2(i,:)), mu_E, norm(r2(i,:)), []);
end

hold on
plot(T,specificEnergy1);
plot(T,specificEnergy2);
hold off
xlabel('Time [s]'); ylabel('$\epsilon [km^2/s^2]$');
title('Specific Energy');
legend('data perturbed', 'ideal data');
grid on;

figure(3)
h_v = zeros(length(T), 3);
h_module = zeros(length(T), 1);
h_v2 = zeros(length(T2), 3);
h_module2 = zeros(length(T2), 1);

for i = 1:length(T)
    [h_v(i,:), h_module(i), ~] = angular_momentum(r(i,:), v(i,:));
end
for i = 1:length(T2)
    [h_v2(i,:), h_module2(i), ~] = angular_momentum(r2(i,:), v2(i,:));
end
subplot(2,1,1)
hold on
plot(T, h_module, 'b-', 'LineWidth', 1.5)
plot(T2, h_module2, 'r--', 'LineWidth', 1.5)
hold off
ylabel('Angular Momentum Magnitude [km^2/s]')
legend('Perturbed', 'Ideal')
grid on

subplot(2,1,2)
hold on
plot(T, h_v(:,3)./h_module, 'b-', 'LineWidth', 1.5) % Normalized z-component
plot(T2, h_v2(:,3)./h_module2, 'r--', 'LineWidth', 1.5)
hold off
xlabel('Time [s]')
ylabel('h_z/h [-]')
legend('Perturbed', 'Ideal')
grid on

figure(4)
% Perturbed data (existing)
[e_v, e, e_unit] = eccentricity(r,v,h_v,mu_E);
e_module = ones(length(T)).*e;

% Ideal data (new)
[e_v2, e2, e_unit2] = eccentricity(r2,v2,h_v2,mu_E);
e_module2 = ones(length(T2)).*e2;

hold on;
% Perturbed
plot(T,e_v(:,1), 'b-');
plot(T,e_v(:,2), 'r-');
plot(T,e_v(:,3), 'g-');
plot(T,e_module, 'k-', 'LineWidth', 2);

% Ideal  
plot(T2,e_v2(:,1), 'b--');
plot(T2,e_v2(:,2), 'r--');
plot(T2,e_v2(:,3), 'g--');
plot(T2,e_module2, 'k--', 'LineWidth', 2);

hold off;
xlabel('Time [s]'); ylabel('e [-]');
title('Eccentricity Vector');
legend('x_p','y_p','z_p','module_p', 'x','y','z','module');
grid on;

figure(5)
% Perturbed data (existing)
v_radial = zeros(length(T), 1);
v_transversal = zeros(length(T), 1);
for i = 1:length(T)
    [radial_vector, ~, transversal_vector, ~] = getvelcomponents(h_v(i,:), v(i,:), r(i,:));
    v_radial(i) = norm(radial_vector);
    v_transversal(i) = norm(transversal_vector);
end

% Ideal data (new)
v_radial2 = zeros(length(T2), 1);
v_transversal2 = zeros(length(T2), 1);
for i = 1:length(T2)
    [radial_vector2, ~, transversal_vector2, ~] = getvelcomponents(h_v2(i,:), v2(i,:), r2(i,:));
    v_radial2(i) = norm(radial_vector2);
    v_transversal2(i) = norm(transversal_vector2);
end

hold on;
% Perturbed
plot(T, v_transversal, 'b-', 'LineWidth', 1.5);
plot(T, v_radial, 'r-', 'LineWidth', 1.5);

% Ideal
plot(T2, v_transversal2, 'b--', 'LineWidth', 1.5);
plot(T2, v_radial2, 'r--', 'LineWidth', 1.5);

hold off;
xlabel('Time [s]'); ylabel('v [km/s]');  % Fixed unit from m/s to km/s
title('Velocity components');
legend('Transversal_p', 'Radial_p', 'Transversal', 'Radial');
grid on;

figure(6)
% Perturbed data (existing)
dotproduct = zeros(length(T), 1);
for i = 1:length(T)
    dotproduct(i) = dot(h_v(i,:), e_v(i,:));
end

% Ideal data (new)
dotproduct2 = zeros(length(T2), 1);
for i = 1:length(T2)
    dotproduct2(i) = dot(h_v2(i,:), e_v2(i,:));
end

hold on;
plot(T, dotproduct, 'b-', 'LineWidth', 1.5);
plot(T2, dotproduct2, 'r--', 'LineWidth', 1.5);
hold off;
xlabel('Time [s]'); ylabel('[-]');
title('e \cdot h dot product');
legend('Perturbed', 'Ideal');
grid on;

% Add a check for orthogonality (should be near zero for ideal case)
fprintf('Mean e-h dot product (Perturbed): %e\n', mean(abs(dotproduct)));
fprintf('Mean e-h dot product (Ideal): %e\n', mean(abs(dotproduct2)));

% Add this to verify energy behavior:
figure(7)
energy_diff = specificEnergy1 - mean(specificEnergy1);
energy_diff2 = specificEnergy2 - mean(specificEnergy2);
plot(T, energy_diff, T2, energy_diff2)
ylabel('Energy Deviation [km^2/s^2]')
xlabel('Time [s]')
title('Energy Conservation')
legend('Perturbed', 'Ideal')
grid on
