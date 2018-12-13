function W = calc_W( L , dataPoints , params)
%calc_W: Calculates the coefficient matrix W that satisfies xdot = W*psi(x,u)
%   Detailed explanation goes here

n = params.n;       % dimension of state, x
nzeta = params.nzeta;   % dimension of state w/ delays, zeta
N = params.Np;       % length of the basis
K = params.K; % size(dataPoints,1);     % total number of datapoints

Ldiag = kron( ones(K,1) , L');    % diagonally stack the transpose of L

% evaluate the basis jacobian at each point and stack the result
dpsi_dx = zeros(K*N, n);
for i = 1 : K
    zeta = dataPoints( i , 1:nzeta )';
%     u = dataPoints( i , (nzeta+1):end )';
    dpsi_dx( (i-1)*N+1 : i*N , : ) =  [ jacobianLift(zeta) ; zeros(params.p , params.n)];
end

W = dpsi_dx \ Ldiag;

end

