function required_a = RepeatingGroundTrack(m,k,omega,mu)
    % RepeatingGroundTrack calculates the required semi-major axis for a 
    % repeating ground track orbit. The function takes into account the 
    % parameters m and k, which are the main parameters of the ground track, 
    % omega, the angular velocity of the planet in radians per second, and 
    % mu, the gravitational parameter of the planet.
    %
    % PROTOTYPE
    %   required_a = RepeatingGroundTrack(m, k, omega, mu)
    %
    % Inputs:
    %   m - integer representing the number of repeat cycles
    %   k - integer representing the number of ground track repeats
    %   omega - angular velocity of the planet (rad/s)
    %   mu - gravitational parameter of the planet (km^3/s^2)
    %
    % Outputs:
    %   required_a - the calculated semi-major axis required for the 
    %   repeating ground track (meters)
    %
    % CONTRIBUTORS:
    %   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
    % 
    % VERSIONS
    %   2025-11: First version
    % 
    % -------------------------------------------------------------------------

    n = omega * (k/ m);
    
    %T = (2*pi/omega)*m/k;
    required_a = (mu /n^2)^(1/3);
end
