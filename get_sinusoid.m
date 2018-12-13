function [ sinusoid ] = get_sinusoid( x, multiplier )
%get_sinusoid: builds a sinusoid from symbolic vector x and a vector of multipliers
%   e.g. x = []

n = length(multiplier); % vector of multipliers

sinusoid = sym(1);  % initialize as a symbolic variable
for i = 1 : n/2
    if multiplier(i) ~= 0
        sinusoid = sinusoid * sin(2*pi*multiplier(i)*x(i));
    end
end
for j = n/2 + 1 : n
    if multiplier(j) ~= 0
        sinusoid = sinusoid * cos(2*pi*multiplier(j)*x(j - n/2));
    end
end


end

