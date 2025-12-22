function PlotPlanet(path,radius)
% PlotPlanet is a function that alows to add a planet or a celestial body
% to a 3D plot. Hold on must be added to the main script.
%
% PROTOTYPE
%   PlotPlanet(path,radius)
%
% INPUT:
%   path[-]          path of the image of the surface to apply [-]
%   radius[1]        radius of the celestial body [km]
%
% CONTRIBUTORS:
%   Federico Masiero
%
% VERSIONS
%   2025-10-09: First version
%
% -------------------------------------------------------------------------

% 1. Uploading and preparing the image
earth = imread(path);
earth = flipud(earth);

% 2. Creation of the sphere
% sphere(n) generates a sphere with (n+1)x(n+1) points
n = 1000;
[X, Y, Z] = sphere(n);

% 3. Create 3D surface adding the texture
surf(X*radius, Y*radius, Z*radius, ...
    'FaceColor', 'texturemap', ...
    'CData', earth, ... 
    'EdgeColor', 'none');

% 4. Perspective
view(3)         

end
