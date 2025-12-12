function [xmin, fval] = find_best_transfer_mod(planet1, planet2, Departure_time_earliest, Departure_time_latest, Arrival_time_earliest, Arrival_time_latest, type, res_t1, res_t2, name)

% Planetary constant of the Sun
mu_Sun = astroConstants(4);

% Planets involved
Earth = planet1;
planet = planet2;

% Conversion to mjd
Departure_time_earliest = date2mjd2000(Departure_time_earliest);
Departure_time_latest   = date2mjd2000(Departure_time_latest);
Arrival_time_earliest   = date2mjd2000(Arrival_time_earliest);
Arrival_time_latest     = date2mjd2000(Arrival_time_latest);

t1 = linspace(Departure_time_earliest, Departure_time_latest, res_t1);
t2 = linspace(Arrival_time_earliest, Arrival_time_latest, res_t2);

Delta_V_tot = zeros(res_t1,res_t2);
ToF = zeros(res_t1,res_t2);

parfor i = 1:res_t1

    for j = 1:res_t2

        [Delta_V_tot(j,i), Dv_1] = Lambert_direct_transfer(Earth,t1(i),planet,t2(j),mu_Sun);

        if norm(Dv_1) > 7
            Delta_V_tot(j,i) = 1000;
        end
        
        ToF(j,i) = t2(j) - t1(i);

    end

    fprintf('%.2f %%\n', (1-i/res_t1+1/res_t1)*100)

end



% Refinement of the solution
[DV_col, j] = min(Delta_V_tot);
[fval, i] = min(DV_col);
j = j(i);

xmin = [t1(i); t2(j)];


if type == 1

    options = optimoptions('fminunc', 'Display', 'none', 'StepTolerance', 1e-16);

    f = @(x) Lambert_transfer_for_fminunc(Earth,x(1),planet,x(2),mu_Sun);

    [xmin, fval] = fminunc(f, xmin, options);

elseif type == 2

    options = optimoptions('fmincon', 'Display', 'none', 'StepTolerance', 1e-16);

    f = @(x) Lambert_transfer_for_fminunc(Earth,x(1),planet,x(2),mu_Sun);

    lb = [xmin(1)-0.01; xmin(2)-0.01];
    ub = [xmin(1)+0.01; xmin(2)+0.01];

    [xmin, fval] = fmincon(f, xmin, [], [], [], [], [], [], @mycon, options);           

end


fprintf('Minimum DV: ');
disp(fval)


% Plot
figure('Name',name)
[X, Y] = meshgrid(t1,t2);
hold on
surf(X,Y,Delta_V_tot,'EdgeColor', 'none')
plot3(xmin(1),xmin(2),fval,'ok','MarkerFaceColor','r','MarkerSize', 5)
view(3)

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);


figure('Name',name)
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
clabel(C, h);
grid on
hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);



figure('Name',name)
[C, h] = contour(X,Y,Delta_V_tot,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
grid on

cb = colorbar;
cb.Label.String = 'ΔV [km/s]';
clim([5 10]);

hold on
plot(xmin(1),xmin(2),'ok','MarkerFaceColor','b','MarkerSize', 5)
 
[C2, h2] = contour(X,Y,ToF,[60, 120, 180, 240, 300],'LineColor', 'k');
clabel(C2, h2);

xlabel('Departure date')
xticks([1185.5 1215.5 1246.5 1276.5 1307.5]);                
xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'});
xtickangle(45);

ylabel('Arrival date')
yticks([1368.5 1399.5 1429.5 1460.5 1491.5 1520.5]);                
yticklabels({'2003 OCt 01','2003 Nov 01','2003 Dec 01','2004 Gen 01','2004 Feb 01','2004 Mar 01'});
ytickangle(45);

% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kep_1, ~] = uplanet (xmin(1), Earth);
[kep_2, ~] = uplanet (xmin(2), planet);

[r1,v1] = kep2car(kep_1, mu_Sun);
[r2,v2] = kep2car(kep_2, mu_Sun);

Delta_t = (xmin(2) - xmin(1)) * 24 * 3600;

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, Delta_t, mu_Sun, 0, 0, 0 );
kep_t = car2kep(r1,v_t1',mu_Sun);

T_Earth = 365.26 * 24 * 3600;
T_planet = 2*pi * sqrt( kep_2(1)^3 / mu_Sun);
T_transf = Delta_t;
T_unused = 2*pi *sqrt( kep_t(1)^3 / mu_Sun) - T_transf;

t0 = 0;

tspan_Earth = linspace(t0,T_Earth,20000);
tspan_planet  = linspace(t0,T_planet,20000);
tspan_transf = linspace(t0,T_transf,20000);
tspan_unused = linspace(t0,T_unused,20000);

S_Earth = [r1; v1];
S_planet  = [r2; v2];
S_transf = [r1; v_t1'];
S_unused = [r2; v_t2'];


% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Propagation of Earth orbit
[ ~ , S_Earth ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_Earth, S_Earth, options );

% Propagation of Earth orbit
[ ~ , S_planet ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_planet, S_planet, options );

% Propagation of Earth orbit
[ ~ , S_transf ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_transf, S_transf, options );

% Propagation of Earth orbit
[ ~ , S_unused ] = ode113( @(t,y) ode_2bp( t, y, mu_Sun ), tspan_unused, S_unused, options );


R_T = 15000000;
[X_T, Y_T, Z_T] = sphere(50);
X_T = R_T * X_T;
Y_T = R_T * Y_T;
Z_T = R_T * Z_T;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

R_M = 15000000;
[X_M, Y_M, Z_M] = sphere(50);
X_M = R_M * X_M;
Y_M = R_M * Y_M;
Z_M = R_M * Z_M;
% mars = imread('MarsTexture.jpg');
% mars = flipud(mars);

R_S = 12000000;
[X_S, Y_S, Z_S] = sphere(50);
X_S = R_S * X_S;
Y_S = R_S * Y_S;
Z_S = R_S * Z_S;
% sun = imread('SunTexture.jpg');
% sun = flipud(sun);


figure('Name',name)

grid on
hold on
axis equal

surf(X_T + r1(1), Y_T + r1(2), Z_T + r1(3), 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
surf(X_M + r2(1), Y_M + r2(2), Z_M + r2(3), 'FaceColor', 'r', 'EdgeColor', 'none');
surf(X_S        , Y_S        , Z_S        , 'FaceColor', 'y', 'EdgeColor', 'y');

plot3(S_Earth(:,1),S_Earth(:,2),S_Earth(:,3),'-b','LineWidth',3)
plot3(S_planet(:,1) ,S_planet(:,2), S_planet(:,3),'-r','LineWidth',3)
plot3(S_transf(:,1),S_transf(:,2),S_transf(:,3),'-g','LineWidth',3)

xlabel('x [km]')
ylabel('y [km]')
view(3)


end














%%

function dy = ode_perturbed_2bp( ~, y, mu, R, J2 )
%ode_perturbed_2bp ODE system for the two-body problem (Keplerian motion)
%considering also the perturbation of Earth's gravity field
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% R[1] Radius [L]
% J2[1] second zonal harmonic [-]

%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez, Federico Masiero
%
% VERSIONS
% 2025-10-06: First version
%
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

% Compute the perturbation as an acceleration

a_J2_x = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(1) / rnorm * (5 * r(3)^2 / rnorm^2 - 1));
a_J2_y = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(2) / rnorm * (5 * r(3)^2 / rnorm^2 - 1));
a_J2_z = (3/2) * J2 * mu * R^2 / rnorm^4 * (r(3) / rnorm * (5 * r(3)^2 / rnorm^2 - 3));

a_J2 = [a_J2_x; a_J2_y; a_J2_z];

% Set the derivatives of the state
dy = [ v; 
    (-mu/rnorm^3)*r + a_J2];
end



function dy = ode_2bp( ~, y, mu )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state
dy = [ v
    (-mu/rnorm^3)*r ];
end




function [DV_tot, DV_1, DV_2] = Lambert_direct_transfer(body1,t1,body2,t2,mu)

ToF = (t2 - t1) * 24 * 3600;

[kep_1, ~] = uplanet (t1, body1);
[kep_2, ~] = uplanet (t2, body2);

[r1,v1] = kep2car(kep_1, mu);
[r2,v2] = kep2car(kep_2, mu);

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, ToF, mu, 0, 0, 0 );

DV_1 = v_t1' - v1;
DV_2 = v2 - v_t2';

DV_tot = norm(DV_1) + norm(DV_2);

end


function [c, ceq] = mycon(x)

    % Associa le variabili
    t1 = x(1);
    t2 = x(2);

    % Costanti note
    body1 = 1;   % inserisci
    body2 = 2;   % inserisci
    mu    = astroConstants(4);   % inserisci

    % Richiama la tua funzione
    c = Lambert_transfer_for_fminunc(body1, t1, body2, t2, mu) - 7;
    ceq = [];
end

function DV_tot = Lambert_transfer_for_fminunc(body1,t1,body2,t2,mu)

ToF = (t2 - t1) * 24 * 3600;

[kep_1, ~] = uplanet (t1, body1);
[kep_2, ~] = uplanet (t2, body2);

[r1,v1] = kep2car(kep_1, mu);
[r2,v2] = kep2car(kep_2, mu);

[~,~,~,~,v_t1,v_t2,~,~] = lambertMR( r1, r2, ToF, mu, 0, 0, 0 );

DV_1 = v_t1' - v1;
DV_2 = v2 - v_t2';

DV_tot = norm(DV_1) + norm(DV_2);

end


