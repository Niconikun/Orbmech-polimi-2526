function E = KeplerEquationSolver(t, e, a, mu, t0, E0) % Not complete!!!
% KeplerEquationSolver Solver of the Kepler Equation
% 
% Function that resolve the Kepler Equation and compute E
% 
% PROTOTYPE
%  E = KeplerEquationSolver(t, e, a, mu, t0, E0)
% 
% INPUT
%   t [-]            Time vector
%   e [1]            Eccentricity vector
%   a [1]            Semi major axis
%   mu [1]           Gravitational parameter of the Earth
%   t0 [1]           Initial time
%   E0 [1]           Initial eccentric anomaly
% 
% OUTPUT
%   E [-]            Vector of eccentric anomaly
% 
% CONTRIBUTORS:
%   Muscas Alice
% 
% VERSIONS
%   2025-10-22: First version


% 1. Compute n
n = sqrt(mu / (a^3));

% 2. Compute E0 and t0 if necessary
if nargin < 6
    E0 = n*t0 + (e* sin(n*t0)) / ( 1 - sin(n*t0 + e) + sin(n*t0) );
    if nargin < 5
        t0 = 0;
    end
end

% 3. Compute E
E = zeros(N, 1);
for i = 1:N
    fun = @(E) E - e*sin(E) - E0 + e*sin(E0) - n * (t(i) - t0);
    E(i) = fzero(fun, E0); 
end

end