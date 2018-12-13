function Uvec = solve_KoopmanQP( Px, Py , params )
%solve_KoopmanQP: Summary of this function goes here
%   Detailed explanation goes here

% na = size( A, 1 );
% nx = size( A, 2 );
% nb = size( b, 1 );

nx = params.N^2;

%% Ridge Method

% Istacked = repmat( speye(params.N) , (size( A,1 ) / params.N) , 1);
% 
% tildeA = [ A , -Istacked, zeros(na,1) ;...
%           -A , -Istacked, zeros(na,1) ;...
%           ones(1,params.N^2) , zeros(1,params.N), -1 ;...
%           -ones(1,params.N^2) , zeros(1,params.N), -1];
% tildeb = [b ; -b; 0 ; 0];
% 
% f = zeros( nx + params.N + 1, 1 );
% f( (nx + 1):(nx + params.N) ) = 1;  % penalty error (was 1)
% f( end ) = 0;     % penalty for size of elements in U (was 100)
% 
% [ xout, fval, exitflag ] = linprog(f, tildeA, tildeb);


%% Lasso method
% x is vectorized Koopman operator, decomposed into positive and negative parts of each entry x = [u11+, ..., uNN+, u11-, ... , uNN-]';
% Uvec = M * x, where M subtracts the + and - parts of each entry: uij+ - uij-

M = [speye(params.Np^2) , -speye(params.Np^2)];

PxTPx = Px' * Px;
PxTPy = Px' * Py;
ATA = kron(speye(params.Np) , PxTPx);  % repeat blocks diagonally N times
ATb = reshape(PxTPy, [params.Np^2 , 1]);

% L2 error as cost function
preH = ATA * M;
H = M' * preH;
f = -M' * ATb;

% L1 regularization enforced as constraint
t = params.t;
Aq = [ -speye(2*params.Np^2) ; ones(1 , 2*params.Np^2) ];
bq = [ zeros(2*params.Np^2 , 1) ; t ];

% Solve the quadratic program
[x , results] = quadprog_gurobi( H , f , Aq , bq );       % use gurobi to solve
% options = optimoptions('quadprog', 'Display', 'iter');
% [ x, fval, exitflag ] = quadprog(H, f, Aq, bq, [], [], [], [], [],options);      % use matlab to solve


% Recover Uvec from the optimization variable
xout = M * x;


%% Set output

Uvec = xout;

end

