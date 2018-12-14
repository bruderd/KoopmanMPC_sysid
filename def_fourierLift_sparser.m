function params = def_fourierLift_sparser( params )
%def_fourierLift: Defines the lifting function that lifts state variable x to
% a fourier basis

[n, p, nzeta, maxDegree] = deal(params.n, params.p, params.nzeta, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
xd = sym('xd', [params.nd * params.n, 1]);   % state delays i.e. for 2 delays: [x_i-1, x_i-2]'
ud = sym('ud', [params.nd * params.p, 1]);   % input delays i.e. for 2 delays: [u_i-1, u_i-2]'
zeta = [x ; xd; ud];    % state variable with delays
u = sym('u', [p, 1]);   % input vector u

% matrix of exponents (N x naug). Each row gives exponents for 1 monomial
multipliers = zeros(1,2*nzeta);
for i = 1:maxDegree
   multipliers = [multipliers; partitions(i, ones(1, 2*nzeta))]; 
end

% Number of basis elements, i.e. dimenstion of p(x)
N = nzeta + size(multipliers , 1);

% create vector of sines and cosines with multipliers
fourierBasis = sym('fourierBasis', [N-nzeta,1]);
for i = 1:N-nzeta
    fourierBasis(i,1) = get_sinusoid(zeta, multipliers(i,:));
end

% put the full state at the beginnig of the basis vector
fourierBasis = [zeta ; fourierBasis];

% create the lifting function: x -> p(x)
liftfunName = [ 'stateLift_' , params.basisID , num2str(params.maxDegree) ];
matlabFunction(fourierBasis, 'File', liftfunName, 'Vars', {zeta});

%% define derivative of lifted state with respect to x

dlift = jacobian(fourierBasis,x);
jacfunName = [ 'jacobianLift_' , params.basisID , num2str(params.maxDegree) ];
matlabFunction(dlift, 'File', jacfunName, 'Vars', {zeta});

%% output variables  
params.Basis = fourierBasis;    % symbolic vector of basis functions, p(x)
params.jacobianBasis = dlift;   % symbolic jacobian of the vector of basis monomials
params.N = N;   % dimension of basis (including the state itself)
params.Np = N + p;  % dimension of lifted state
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xd = xd; % symbolic state delays
params.ud = ud; % symbolic input delays
params.liftHandle = liftFunName;    % name of the lifting function

end