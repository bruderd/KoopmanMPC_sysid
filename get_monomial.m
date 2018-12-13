function [ monomial ] = get_monomial( x, exponents )
%get_monomial: builds a monomial from symbolic vector x and a vector of
%exponents
%   e.g. x = [x1 x2]; exponents = [1 2]; =>  monomial = x1^1 * x2^2

n = length(x);

monomial = x(1)^exponents(1);
for i = 2:n
    monomial = monomial * x(i)^exponents(i);
end

end