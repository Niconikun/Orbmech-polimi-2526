function outputArg1 = semilatus_rectum(h,mu,a,e)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    arguments
        h {isscalar}
        mu {isscalar}
        a {isscalar}
        e {isscalar}
        end
    if ~isempty(e) && ~isempty(a)
        outputArg1 = a * (1 - e^2); % Calculate the second output argument
    elseif ~isempty(h) && ~isempty(mu)
        outputArg1 = (h^2) / mu; % Calculate the semi-latus rectum
    else
        error('Insufficient input parameters provided.'); % Handle error case
    end
end