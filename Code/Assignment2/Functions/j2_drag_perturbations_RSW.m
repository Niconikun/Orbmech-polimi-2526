function a_per = j2_drag_perturbations_RSW(t, kep, mu, J2, R, omega, c_d, AoverMass)
% J2_DRAG_PERTURBATIONS_RSW Returns J2 and drag perturbations in RSW frame
%
% INPUT:
%   t         - Time (s)
%   kep       - Keplerian elements [6x1]
%   mu        - Gravitational parameter [km^3/s^2]
%   J2        - J2 harmonic coefficient
%   R         - Planetary radius [km]
%   omega     - Planetary rotation vector [3x1] (rad/s)
%   c_d       - Drag coefficient
%   AoverMass - Area-to-mass ratio [m^2/kg]
%
% OUTPUT:
%   a_per     - Perturbing acceleration in RSW frame [3x1] (km/s^2)

    % Convert Keplerian to Cartesian
    [r, v] = kep2car(kep, mu);
    
    % Get perturbations in Cartesian frame
    % Using your existing ode_2bp_j2_drag function logic
    rnorm = norm(r);
    x = r(1); y = r(2); z = r(3);
    
    % J2 acceleration in Cartesian frame
    factor = (3*J2*mu*(R^2))/(2*(rnorm^4));
    a_j2_cart = factor * [x/rnorm * (5*(z/rnorm)^2 - 1);
                         y/rnorm * (5*(z/rnorm)^2 - 1);
                         z/rnorm * (5*(z/rnorm)^2 - 3)];
    
    % Atmospheric drag in Cartesian frame
    altitude = rnorm - R;
    rho = ExponentialaAtmosphericModel(altitude);  % Your atmospheric model
    
    % Relative velocity
    v_rel = v - cross(omega, r);
    a_drag_cart = -0.5 * c_d * AoverMass * rho * norm(v_rel) * v_rel;
    
    % Total perturbation in Cartesian frame
    a_pert_cart = a_j2_cart + a_drag_cart;
    
    % Convert to RSW frame
    a_per = cartesian_to_RSW(r, v, a_pert_cart);
end