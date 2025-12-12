function [vector, module, unit_vector] = eccentricity(position, velocity, angular_momentum, mu)
%Provides you with the eccentricity vector and its module
    vector = (1./mu).*cross(velocity, angular_momentum) - position./norm(position);
    module = norm(vector);
    unit_vector = vector./module;
end