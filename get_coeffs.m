function CF = get_coeffs( vf, params )
%get_coeffs: Get the coefficients of the symbolic expression vf
%corresponding to monomial basis
%   Detailed explanation goes here

CF = zeros(params.n, params.N);

for i = 1:params.n
    [C, T] = coeffs( vf(i,:) );
    for j = 1:params.N
        for k = 1:length(T)
            if T(k) == params.polyBasis(j)
                CF(i,j) = C(k);
                break
            else
                CF(i,j) = 0;
            end
        end
    end
end


end

