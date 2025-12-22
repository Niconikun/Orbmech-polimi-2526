function [r, v] = keplerianToCartesian(a, e, i, Omega, omega, nu)
    % Convert Keplerian elements to Cartesian coordinates
    % Inputs:
    % a     - semi-major axis (km)
    % e     - eccentricity
    % i     - inclination (radians)
    % Omega  - right ascension of ascending node (radians)
    % omega  - argument of periapsis (radians)
    % nu    - true anomaly (radians)
    %
    % Outputs:
    % r     - position vector in Cartesian coordinates (km)
    % v     - velocity vector in Cartesian coordinates (km/s)

    % Gravitational parameter (mu) for Earth (km^3/s^2)
    mu = 398600;

    % Calculate the distance from the central body
    r_pqw = [(a * (1 - e^2)) / (1 + e * cos(nu)) * cos(nu);
             (a * (1 - e^2)) / (1 + e * cos(nu)) * sin(nu);
             0];

    % Calculate the velocity in the perifocal coordinate system
    v_pqw = [sqrt(mu / a) * sin(nu);
             sqrt(mu / a) * (e + cos(nu));
             0];

    % Rotation matrices
    R3_Omega = [cos(Omega), sin(Omega), 0;
                -sin(Omega), cos(Omega), 0;
                0, 0, 1];

    R1_i = [1, 0, 0;
            0, cos(i), sin(i);
            0, -sin(i), cos(i)];

    R3_omega = [cos(omega), sin(omega), 0;
                -sin(omega), cos(omega), 0;
                0, 0, 1];

    % Combined rotation matrix
    R = R3_Omega * R1_i * R3_omega;

    % Convert to geocentric coordinates
    r = R * r_pqw;
    v = R * v_pqw;
end