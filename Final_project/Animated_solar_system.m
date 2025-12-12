clc
clear
close all

AU = astroConstants(2);
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

xlim([-9 9]*AU)
ylim([-9 9]*AU)
zlim([-4 4]*AU)

% --- Disegno del Sole ---
R_S = 50000000;
[X, Y, Z] = sphere(50);
surf(R_S*X, R_S*Y, R_S*Z, 'FaceColor', 'y', 'EdgeColor', 'none');

% --- Lunghezza scia (numero di punti mantenuti) ---
N = 2000; 

% --- Buffer iniziali vuoti per i pianeti ---
Earth_trail    = nan(N,3);
Mars_trail     = nan(N,3);
Jupiter_trail  = nan(N,3);
Asteroid_trail = nan(N,3);

% --- Oggetti grafici: scie ---
pEarth_trail    = plot3(nan, nan, nan, '.','Color','b','MarkerSize',1);
pMars_trail     = plot3(nan, nan, nan, '.','Color','r','MarkerSize',1);
pJupiter_trail  = plot3(nan, nan, nan, '.','Color','m','MarkerSize',1);
pAsteroid_trail = plot3(nan, nan, nan, '.','Color',[0.8 0.8 0.8],'MarkerSize',1);

% --- Oggetti grafici: teste (posizione attuale) ---
pEarth_head    = plot3(nan, nan, nan, 'o','MarkerFaceColor','b','MarkerEdgeColor','w','MarkerSize',4);
pMars_head     = plot3(nan, nan, nan, 'o','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',4);
pJupiter_head  = plot3(nan, nan, nan, 'o','MarkerFaceColor','m','MarkerEdgeColor','w','MarkerSize',7);
pAsteroid_head = plot3(nan, nan, nan, 'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','w','MarkerSize',2);

% ================================
%       LOOP DI ANIMAZIONE
% ================================
for i = 1:res

    time = Time_interval(i);
    standard_time = mjd20002date(time);
    month = ["January","February","March","April","May","June","July","August","September","October","November","December"];
    year = num2str(standard_time(1));
    day = num2str(standard_time(3));

    standard_time = year + " " + month(standard_time(2)) + " " + day;
    
    title(sprintf('%s', standard_time));

    % --- Calcolo posizioni ---
    kep_Earth    = uplanet(time, Earth);
    kep_Mars     = uplanet(time, Mars);
    kep_Jupiter  = uplanet(time, Jupiter);
    kep_asteroid = ephAsteroids(time, asteroid);

    % Converti kepleriane → coordinate (X,Y,Z)
    rE = kep2car(kep_Earth,    mu_Sun);
    rM = kep2car(kep_Mars,     mu_Sun);
    rJ = kep2car(kep_Jupiter,  mu_Sun);
    rA = kep2car(kep_asteroid, mu_Sun);

    % --- Shift buffers (scie a lunghezza costante N) ---
    Earth_trail    = [Earth_trail(2:end,:);    rE(:)'];
    Mars_trail     = [Mars_trail(2:end,:);     rM(:)'];
    Jupiter_trail  = [Jupiter_trail(2:end,:);  rJ(:)'];
    Asteroid_trail = [Asteroid_trail(2:end,:); rA(:)'];

    % --- Aggiorna scie ---
    pEarth_trail.XData    = Earth_trail(:,1); pEarth_trail.YData = Earth_trail(:,2); pEarth_trail.ZData = Earth_trail(:,3);
    pMars_trail.XData     = Mars_trail(:,1);  pMars_trail.YData  = Mars_trail(:,2);  pMars_trail.ZData  = Mars_trail(:,3);
    pJupiter_trail.XData  = Jupiter_trail(:,1);pJupiter_trail.YData= Jupiter_trail(:,2);pJupiter_trail.ZData= Jupiter_trail(:,3);
    pAsteroid_trail.XData = Asteroid_trail(:,1);pAsteroid_trail.YData=Asteroid_trail(:,2);pAsteroid_trail.ZData=Asteroid_trail(:,3);

    % --- Aggiorna teste ---
    pEarth_head.XData = rE(1); pEarth_head.YData = rE(2); pEarth_head.ZData = rE(3);
    pMars_head.XData  = rM(1); pMars_head.YData  = rM(2); pMars_head.ZData  = rM(3);
    pJupiter_head.XData = rJ(1); pJupiter_head.YData = rJ(2); pJupiter_head.ZData = rJ(3);
    pAsteroid_head.XData = rA(1); pAsteroid_head.YData = rA(2); pAsteroid_head.ZData = rA(3);

    pause(0.01);

    drawnow limitrate nocallbacks
end
