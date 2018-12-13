function Ainv = dinv( A, damping, params )
%dinv: Returns the pseudoinverse of A
%   Damps A to make it better conditioned for inversion.
%   damping :   the amount to damp the matrix (usually around 1e-3)

[U,S,V] = svds(A, params.N);

for i = 1 : min(size(S,2), size(S,1))
   S(i,i) = max(abs(S(i,i)), damping);
end

Sinv = pinv(S);

Ainv = V * Sinv * U';


end

