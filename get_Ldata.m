function Ldata = get_Ldata( U, params )
%get_Ldata: Calculate the approximate infinitesimal generator for the
%Koopman operator U
%   U = expm(L Ts). This function solves for L

% % eliminate complex eigenvalues of U so that it can be inverted
% [V, D, W] = eig(U);
% Dreal = D;
% for i = 1 : size(Dreal,2)
%    if Dreal(i,i) == 0
%        Dreal(i,i) = 1e-8;
%    end
% end
% logU = V * logm(Dreal) * W';
% 
% Ldata = (1/params.Ts) * logU;
% Ldata = real(Ldata);    % ignore complex components


% % Try using the SVD trick to take the matrix logarithm
% svals = svd(U);
% nzsvals = (svals > 1e-6);   % round close ones down to zero
% numnz = nnz(nzsvals); % number of nonzero singular values
% [u,s,v] = svds(U, numnz);
% logU = v * logm(s) * u';
% Ldata = (1/params.Ts) * logU;

Ldata = (1/params.Ts) * logm(U);

end