%% Exercise 1a

clc
clear
close all


% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_Earth = astroConstants(13);
mu_Sun = astroConstants(4);
AU = astroConstants(2);

v_inf = [15.1; 0; 0];
r = [1; 0; 0] * AU;

Delta = 9200;

[a, delta, e, rp] = solve_hyperbola_2d(norm(v_inf),Delta,mu_Earth);


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In front of the planet
u = [0; 0; -1];
v_front = rodrigues_rotation(u,v_inf,delta);

% Behind of the planet
u = [0; 0; 1];
v_behind = rodrigues_rotation(u,v_inf,delta);

% Under of the planet
u = [0; -1; 0];
v_under = rodrigues_rotation(u,v_inf,delta);


% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_T = sqrt(mu_Sun / AU) * [0; 1; 0];

V_IN = V_T + v_inf;

% In front of the planet
V_FRONT = V_T + v_front;

% Behind of the planet
V_BEHIND = V_T + v_behind;

% Under of the planet
V_UNDER = V_T + v_under;


% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% Exercise 1b

clc
clear
close all


mu_Earth = astroConstants(13);
mu_Sun = astroConstants(4);
AU = astroConstants(2);

v_inf = [15.1; 0; 0];
r_Earth = [1; 0; 0] * AU;

% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I choose the incoming asymptote in front of the Earth


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining different values of impact parameter
Delta = linspace(9200, 13200, 5);

% Solving different hyperbolas
[a, delta, e, rp] = arrayfun(@(x) solve_hyperbola_2d(norm(v_inf), x, mu_Earth), Delta);

% Earth's velocity around Sun
V_T = sqrt(mu_Sun / AU) * [0; 1; 0];

% Rotation axis of v_inf_meno
u = [0; 0; -1];

% v_inf using rotation of the vector and Rodrigues method
v_front = zeros(3, length(delta));
for i = 1:length(delta)
    v_front(:,i) = rodrigues_rotation(u, v_inf, delta(i));
end

% Orbital parameter of the hyperbolas
kep = zeros(length(delta),6);

for i = 1:length(delta)
    [i_front, Omega_front, omega_front] = plane(v_inf,v_front(:,i),delta(i));
    kep(i,:) = [a(i), e(i),i_front, Omega_front, omega_front, 0];
end


% Earth image
R_T = astroConstants(23);
[X, Y, Z] = sphere(50);
X_T = R_T*X;
Y_T = R_T*Y;
Z_T = R_T*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

% Plot of hyperbolas
figure('Name','Hyperbolas')
hold on

%Trajectories
for i = 1:length(Delta)
    my_plotOrbit(kep(i,:), mu_Earth, 1000, -1.6, +1.6, '')
end

%Percientres
for i = 1:length(Delta)
    my_plotOrbit(kep(i,:), mu_Earth, 1, 0, 0, 'o')
end

% Earth's velocity
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)

% Earth image
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');

% Settings
legend('\Delta = 9200 [km]', '\Delta = 10200 [km]', '\Delta = 11200 [km]', '\Delta = 12200 [km]', '\Delta = 13200 [km]', 'Interpreter', 'Tex')
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
xlim([-4 4] * 1e4)
ylim([-3 3] * 1e4)


% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_IN = V_T + v_inf;

KEP_IN = car2kep_mod(r_Earth,V_IN,mu_Sun);


% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_OUT = V_T + v_front;

KEP_OUT = zeros(length(Delta), 6);

for i = 1:length(Delta)
    KEP_OUT(i,:) = car2kep_mod(r_Earth,V_OUT(:,i),mu_Sun);
end


% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_S = 10000000;

X_S = R_S*X;
Y_S = R_S*Y;
Z_S = R_S*Z;


% Plot in front of the planet 
figure('Name','Heliocentric trajectory')
hold on

% Outcoming trajectories
for i = 1:length(Delta)
    my_plotOrbit(KEP_OUT(i,:), mu_Sun, 1000, KEP_OUT(i,6), pi + KEP_OUT(i,6), '')
end

% Incoming trajectory
my_plotOrbit(KEP_IN, mu_Sun, 1000, KEP_IN(6)-3/2*pi, KEP_IN(6), '')

% Fly-by point
my_plotOrbit(KEP_IN, mu_Sun, 1, KEP_IN(6), KEP_IN(6), 'ok')

% Sun
surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');

% Settings
legend('\Delta = 9200 [km]', '\Delta = 10200 [km]', '\Delta = 11200 [km]', '\Delta = 12200 [km]', '\Delta = 13200 [km]', 'Interpreter', 'Tex')
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


% Defining different values of impact parameter
Delta_plot = linspace(1, 39, 100);

% Solving different hyperbolas
[~, delta_plot, ~, rp_plot] = arrayfun(@(x) solve_hyperbola_2d(norm(v_inf), x, mu_Earth), Delta_plot*R_T);


% Plot of the trends
figure('Name','Trends of altitude and turning angle')

yyaxis left
plot(Delta_plot,rp_plot/R_T - 1,'b','LineWidth',2)
ylabel('Fly by minimum altitude')
ax.YColor = 'b';
ylim([0 40]); 

yyaxis right
plot(Delta_plot,delta_plot * 180/pi,'r','LineWidth',2)
ylabel('Turning angle \delta', 'Interpreter','Tex')
ax.YColor = 'r';
ylim([0 35]); 

xlabel('Impact parameter \Delta', 'Interpreter','Tex')
grid on
legend('Flyby minimum altitude','Turning angle')


%% Exercise 2

clc
clear
close all

mu_Earth = astroConstants(13);
mu_Sun = astroConstants(4);
AU = astroConstants(2);

V_MINUS = [31.5; 5.2; 0.0];
V_PLUS = [36.0; 0.0; 0.0];
R_EARTH = [0; -1; 0] * AU;

% Earth's velocity around Sun
V_T = sqrt(mu_Sun / AU) * [0; 1; 0];


% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_EARTH = sqrt(mu_Sun / AU) * [1; 0; 0];

v_inf_minus = V_MINUS - V_EARTH;

v_inf_plus = V_PLUS - V_EARTH;


% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = acos( dot(v_inf_plus,v_inf_minus) / ( norm(v_inf_plus)*norm(v_inf_minus) ) );


% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_T = astroConstants(23);

v_in_norm = norm(v_inf_minus);
v_out_norm = norm(v_inf_plus);

% Searching radius of pericentre
fun = @(r_p) asin( 1 / ( 1 + r_p*(v_out_norm^2)/mu_Earth ) ) + asin( 1 / ( 1 + r_p*(v_in_norm^2)/mu_Earth ) ) - delta;
x0 = R_T + 1000;
options = optimset('TolX', 1e-12, 'TolFun', 1e-16, 'MaxIter', 10000, 'MaxFunEvals', 50000 ); 
r_p = fzero(fun,x0,options);

if r_p > R_T
    disp('The radius of the percienter is grater than the radius of the Earth')
else
    disp('The radius of the percienter is lower than the radius of the Earth, so the fly-by is not feasible')
end


% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_p_out = sqrt(v_out_norm^2 + 2*mu_Earth/r_p);
v_p_in  = sqrt(v_in_norm^2  + 2*mu_Earth/r_p);

dv_p = abs(v_p_out - v_p_in);

dv_fly_by = norm(V_PLUS - V_MINUS);


% Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_in  = - mu_Earth / v_in_norm^2;
a_out = - mu_Earth / v_out_norm^2;

e_in  = 1 - r_p / a_in;
e_out = 1 - r_p / a_out;

delta_in  = 2 * asin( 1 / e_in );
delta_out = 2 * asin( 1 / e_out );

[i_in, Omega_in, omega_in]    = plane(v_inf_minus,v_inf_plus,delta_out);
[i_out, Omega_out, omega_out] = plane(v_inf_minus,v_inf_plus,delta_out);

kep_in  = [a_in,  e_in,  i_in,  Omega_in,  omega_in,  0];
kep_out = [a_out, e_out, i_out, Omega_out, omega_out, 0];

KEP_BEFORE = car2kep_mod(R_EARTH,V_MINUS,mu_Sun);
KEP_AFTER  = car2kep_mod(R_EARTH,V_PLUS ,mu_Sun);

c_in  = abs( a_in  * e_in );
c_out = abs( a_out * e_out );

h_dir = cross(v_inf_minus,v_inf_plus) / norm( cross(v_inf_minus,v_inf_plus) );

r_p_dir = -rodrigues_rotation(h_dir,v_inf_plus,pi/2-delta_out/2) / v_out_norm;

r_p_in  = c_in  * r_p_dir;
r_p_out = c_out * r_p_dir;

asymptote_in  = r_p_in  + v_inf_minus * [-1 -10000];
asymptote_out = r_p_out + v_inf_plus  * [ 1  10000];

tratteggio = r_p_in * [1 -2.5];

% Plot planet 
[X, Y, Z] = sphere(50);
X_T = R_T*X;
Y_T = R_T*Y;
Z_T = R_T*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);

R_S = 10000000;

X_S = R_S*X;
Y_S = R_S*Y;
Z_S = R_S*Z;

figure('Name','Power gravity assist')

subplot(1,2,1)
title('Heliocentric frame')
hold on
my_plotOrbit(KEP_BEFORE, mu_Sun, 1000, KEP_BEFORE(6)-0.7*pi, KEP_BEFORE(6)          , 'b')
my_plotOrbit(KEP_AFTER , mu_Sun, 1000, KEP_AFTER(6)        , KEP_AFTER(6) + 1.5*pi/2, 'r')
my_plotOrbit(KEP_BEFORE, mu_Sun, 1   , KEP_BEFORE(6)       , KEP_BEFORE(6)          , 'ok')
surf(X_S, Y_S, Z_S, 'FaceColor', 'y', 'EdgeColor', 'none');
axis equal
grid on
legend('Before flyby', 'After flyby')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(3)

subplot(1,2,2)
title('Planetocentric frame')
hold on
my_plotOrbit(kep_in , mu_Earth, 1000, -1.9,    0, 'b')
my_plotOrbit(kep_out, mu_Earth, 1000, 0   , +1.9, 'r')
plot3(asymptote_in(1,:),asymptote_in(2,:),asymptote_in(3,:),'LineWidth',1.5)
plot3(asymptote_out(1,:),asymptote_out(2,:),asymptote_out(3,:),'LineWidth',1.5)
plot3(tratteggio(1,:),tratteggio(2,:),tratteggio(3,:),'--k','LineWidth',1.5)
my_plotOrbit(kep_out, mu_Earth, 1   , 0   ,    0, 'ob')
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
axis equal
grid on
legend('Incoming hyperbola', 'Outcoming hyperbola','Asymptote of incoming hyperbola','Asymptote of outcoming hyperbola')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(3)


%% Additional exercise - part 1

clc
clear
close all

AU       = astroConstants(2);
mu_Sun   = astroConstants(4);
mu_Earth = astroConstants(13);

delta = 25.8230 * pi/180;
v_inf_plus = [10.0856; -6.4936; 0];
R_T = [ 0; -1; 0] * AU;


% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_inf_norm = norm(v_inf_plus);

a = - mu_Earth / v_inf_norm^2;

Delta = - a / tan( delta/2 );


% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = [0; 0; 1];

v_inf_minus = rodrigues_rotation(u, v_inf_plus, delta);


% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_T = sqrt( mu_Sun / AU ) * [1; 0; 0];

V_MINUS = V_T + v_inf_minus;


% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = 1 / sin( delta/2 );

[i, Omega, omega] = plane(v_inf_minus,v_inf_plus,delta);

kep = [ a, e, i, Omega, omega, 0];


r_t = astroConstants(23);

[X, Y, Z] = sphere(50);
X_T = r_t*X;
Y_T = r_t*Y;
Z_T = r_t*Z;
earth = imread('EarthTexture.jpg');
earth = flipud(earth);


figure('Name','Additional exercise - Part 1')

title('Planetocentric frame')
hold on
surf(X_T, Y_T, Z_T, 'FaceColor', 'texturemap', 'CData', earth, 'EdgeColor', 'none');
my_plotOrbit(kep, mu_Earth, 1000, -1.5, +1.5, 'b')
my_plotOrbit(kep, mu_Earth, 1, 0, 0, 'ob')
quiver3(0,0,0,V_T(1)*1000,V_T(2)*1000,V_T(3)*1000,'k','LineWidth',2)
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


%% Additional exercise - Part 2

clc
clear
close all

AU = astroConstants(2);
mu_Sun = astroConstants(4);
mu_Earth = astroConstants(13);

v_inf_plus = [10.0856; -6.4936; 0];
v_inf_minus = [5.6934; 0; 0];
R_T = [-1; 0; 0] * AU;


% Question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_T = sqrt( mu_Sun / AU ) * [0; -1; 0];

V_MINUS = V_T + v_inf_minus;


% Question 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_plus_norm = norm(v_inf_plus);
v_minus_norm = norm(v_inf_minus);

delta = acos( dot(v_inf_plus,v_inf_minus) / (v_plus_norm*v_minus_norm) );


% Question 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_plus = - mu_Earth / v_plus_norm^2;


% Question 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_t = astroConstants(23);

fun = @(r_p) asin( 1 / ( 1 + r_p*(v_plus_norm^2)/mu_Earth ) ) + asin( 1 / ( 1 + r_p*(v_minus_norm^2)/mu_Earth ) ) - delta;
x0 = r_t + 6000;
options = optimset('TolX', 1e-12, 'TolFun', 1e-16, 'MaxIter', 10000, 'MaxFunEvals', 50000 ); 
r_p = fzero(fun,x0,options);

% Using conservation of energy
v_p_minus = sqrt( v_minus_norm^2   +   2 * mu_Earth / r_p );
v_p_plus  = sqrt( v_plus_norm^2    +   2 * mu_Earth / r_p );

Dv = v_p_plus - v_p_minus;


% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rp = linspace(r_t+1000,r_t+20000,1000000);

y = arrayfun(fun,rp);

figure('Name','Function to minimise')
hold on
plot(rp,y,'DisplayName','function')
plot(x0,fun(x0),'*r','DisplayName','initial guess')
plot(r_p,fun(r_p),'or','DisplayName','found value')
grid on
legend('show')
xlabel('Radius of pericenter [km]')
ylabel('Value of the function')


%% Test car2kep

kep_i = [20000, 0.2, pi, 0, pi+pi/3, 0];

[r, v] = kep2car(kep_i,mu_Earth);

r(3) = 0;
v(3) = 0;

car2kep(r,v,mu_Earth)



%% FUNCTIONS

function [a, delta, e, rp] = solve_hyperbola_2d(v_inf_meno,Delta,mu)

a = -mu / v_inf_meno^2;

delta = 2 * atan(-a / Delta );

e = 1 / sin(delta / 2);

rp = - mu / v_inf_meno^2 * ( 1 - e );


end




function v_rot = rodrigues_rotation(u,v,delta)

v_rot = v*cos(delta) + cross(u,v)*sin(delta) + u*dot(u,v)*(1-cos(delta));

end





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




function [r_SOI] = SOI(M,m,distance)

r_SOI = distance*(m/M)^(2/5);

end




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