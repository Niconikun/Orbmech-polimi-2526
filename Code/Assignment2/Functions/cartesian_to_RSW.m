function a_RSW = cartesian_to_RSW(r, v, a_cart)
% CARTESIAN_TO_RSW Converts acceleration from Cartesian to RSW frame
%
% PROTOTYPE
%   a_RSW = cartesian_to_RSW(r, v, a_cart)
%
% INPUT:
%   r      - Position vector in inertial frame [3x1] (km)
%   v      - Velocity vector in inertial frame [3x1] (km/s)
%   a_cart - Acceleration vector in inertial frame [3x1] (km/s^2)
%
% OUTPUT:
%   a_RSW  - Acceleration in RSW frame [3x1] (km/s^2)

% CONTRIBUTORS:
%   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
% 
% VERSIONS
%   2025-11: First version
% 
% -------------------------------------------------------------------------

    % Unit vectors in RSW frame
    r_norm = norm(r);
    R_hat = r / r_norm;  % Radial direction
    
    h_vec = cross(r, v);  % Angular momentum vector
    W_hat = h_vec / norm(h_vec);  % Out-of-plane direction
    
    S_hat = cross(W_hat, R_hat);  % Transversal direction
    
    % Transformation matrix from inertial to RSW
    Q = [R_hat, S_hat, W_hat]';
    
    % Transform acceleration
    a_RSW = Q * a_cart;
end