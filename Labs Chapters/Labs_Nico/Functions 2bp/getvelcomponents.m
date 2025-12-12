function [radial_vector,unit_radial_vector, transversal_vector, unit_transversal_vector] = getvelcomponents(h,v,r)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    h {isvector}
    v {isvector}
    r {isvector}

end

    unit_radial_vector = r./norm(r);

    radial_vector = v.*unit_radial_vector;
    unit_h = h./norm(h);
    unit_transversal_vector = cross(unit_h, unit_radial_vector);
    transversal_vector = v.*unit_transversal_vector;

end