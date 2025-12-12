function [vector, module, unit_vector] = angular_momentum(vector1,vector2)

% Vector and module of angular momentum
module = norm(cross(vector1, vector2));
vector = cross(vector1, vector2);
unit_vector = vector./module;
%   Detailed explanation goes here

end
