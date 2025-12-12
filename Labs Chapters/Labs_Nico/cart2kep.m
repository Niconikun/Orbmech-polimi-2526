function kep = cart2kep(r, v, mu)
% CART2KEP Convert Cartesian state vectors to Keplerian elements
%
% Inputs:
%   r  - Position vector [x; y; z] in km
%   v  - Velocity vector [vx; vy; vz] in km/s
%   mu - Standard gravitational parameter (km^3/s^2). Default: Earth (398600.4418)
%
% Outputs:
%   kep - Keplerian elements [a; e; i; Omega; omega; M] where:
%         a: Semi-major axis (km)
%         e: Eccentricity
%         i: Inclination (rad)
%         Omega: RAAN (rad)
%         omega: Argument of periapsis (rad)
%         M: Mean anomaly (rad)

    if nargin < 3
        mu = 398600.4418; % Earth's gravitational parameter
    end

    r_norm = norm(r);
    v_norm = norm(v);

    % Specific angular momentum
    h = cross(r, v);
    h_norm = norm(h);

    % Eccentricity vector
    e_vec = ((v_norm^2 - mu/r_norm) * r - dot(r, v) * v) / mu;
    e = norm(e_vec);

    % Semi-major axis
    E = v_norm^2/2 - mu/r_norm;
    if e ~= 1
        a = -mu/(2*E);
    else
        a = inf; % Parabolic case
    end

    % Inclination
    i = acos(h(3)/h_norm);

    % RAAN (Ω)
    n = cross([0; 0; 1], h);
    n_norm = norm(n);
    if n_norm ~= 0
        Omega = acos(n(1)/n_norm);
        if n(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0; % Equatorial orbit
    end

    % Argument of periapsis (ω)
    if n_norm ~= 0
        omega = acos(dot(n, e_vec)/(n_norm * e));
        if e_vec(3) < 0
            omega = 2*pi - omega;
        end
    else
        % Equatorial orbit
        omega = atan2(e_vec(2), e_vec(1));
        if omega < 0
            omega = omega + 2*pi;
        end
    end

    % True anomaly (θ)
    if e > 0
        theta = acos(dot(e_vec, r)/(e * r_norm));
        if dot(r, v) < 0
            theta = 2*pi - theta;
        end
    else
        % Circular orbit
        theta = acos(dot(n, r)/(n_norm * r_norm));
        if r(3) < 0
            theta = 2*pi - theta;
        end
    end

    % Mean anomaly (M)
    if e < 1
        % Elliptical
        E = 2 * atan2(sqrt(1-e) * tan(theta/2), sqrt(1+e));
        M = E - e*sin(E);
    elseif e == 1
        % Parabolic
        M = tan(theta/2) + (tan(theta/2))^3/3;
    else
        % Hyperbolic
        F = 2 * atanh(sqrt((e-1)/(e+1)) * tan(theta/2));
        M = e*sinh(F) - F;
    end

    % Ensure M is in [0, 2π)
    M = mod(M, 2*pi);

    kep = [a; e; i; Omega; omega; M];
end