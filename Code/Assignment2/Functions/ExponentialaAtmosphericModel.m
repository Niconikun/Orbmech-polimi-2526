function density = ExponentialaAtmosphericModel(altitude)
    % Exponential Atmospheric Model Documentation
    % 
    % This function calculates the atmospheric density at a given altitude 
    % using an exponential decay model. The model is based on predefined 
    % scale heights and nominal densities at various base altitudes.
    %
    % Inputs:
    %   altitude - The altitude (in meters) at which to calculate the 
    %              atmospheric density.
    %
    % Outputs:
    %   density - The calculated atmospheric density (in kg/m^3) at the 
    %             specified altitude. If the altitude is above the maximum 
    %             defined altitude, the density is returned as 0.
    %
    % The model uses two arrays: 'scale_heights' and 'nominal_density', 
    % which contain the scale heights and corresponding nominal densities 
    % for different altitude ranges. The function iterates through the 
    % 'base_altitude' array to find the appropriate range for the input 
    % altitude and applies the exponential formula to compute the density.
    %
    %
    % CONTRIBUTORS:
    %   Muscas Alice, Masiero Federico, Karthikeyan Prthik Nandhan, Nicolás Sepúlveda
    % 
    % VERSIONS
    %   2025-12-22: First version
    % -------------------------------------------------------------------------


    scale_heights = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 ... 
        5.877 7.263 9.473 12.636 16.149 22.523 29.470 37.105 45.546 53.628 ...
        53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00];
    
    nominal_density = [1.225 3.899e-2 1.774e-2 3.972e-3 1.057e-3 3.206e-4 ...
        8.770e-5 1.905e-5 3.396e-6 5.297e-7 9.661e-8 2.438e-8 8.484e-9 ...
        3.845e-9 2.070e-9 5.464e-10 2.789e-10 7.248e-11 2.418e-11 9.158e-12 ...
        3.725e-12 1.585e-12 6.967e-13 1.454e-13 3.614e-14 1.170e-14 5.245e-15 3.019e-15];

    base_altitude = [0 25 30 40 50 60 70 80 90 10 110 120 130 140 150 180 ...
        200 250 300 350 400 450 500 600 700 800 900 1000];

    for i=1:length(base_altitude)
        if altitude >= base_altitude(i) && altitude < base_altitude(i+1)
            density = nominal_density(i) * exp(-(altitude - base_altitude(i)) / scale_heights(i));
            break;
        elseif altitude > base_altitude(end)
            density = 0;
            break;
        end
    end
end