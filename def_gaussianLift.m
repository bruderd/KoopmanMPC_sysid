function params = def_gaussianLift( params )
%def_thinplateLift: Defines the lifting function that lifts state variable x to
% a set of thin plate spline radial basis functions 

[n, p, nzeta, maxDegree] = deal(params.n, params.p, params.nzeta, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
xd = sym('xd', [params.nd * params.n, 1]);   % state delays i.e. for 2 delays: [x_i-1, x_i-2]'
ud = sym('ud', [params.nd * params.p, 1]);   % input delays i.e. for 2 delays: [u_i-1, u_i-2]'
zeta = [x ; xd; ud];    % state variable with delays
u = sym('u', [p, 1]);   % input vector u

% Number of basis elements, i.e. dimenstion of p(x)
N = nzeta + maxDegree; 

% create basis functions with random centers in interval [-params.scale, params.scale]
psi = sym('psi', [maxDegree , 1]);
zeta0 = params.scale * (2*rand([nzeta,maxDegree]) - 1); % columns are random centers
for i = 1 : maxDegree
   radius = norm( zeta - zeta0(:,i) );
   psi(i,:) = exp(-( 1 * radius )^2) ;
%    % I think this might work faster
%    radius = sum( ( zeta - zeta0(:,i) ).^2 );
%    psi(i,:) = exp( -radius );
end

% define basis vector, putting the full state at the beginnig of the basis vector
Basis = [zeta ; psi];

% create the lifting function: x -> p(x)
liftfunName = [ 'stateLift_' , params.basisID , num2string(params.maxDegree) ];
matlabFunction(Basis, 'File', liftfunName, 'Vars', {zeta});

%% define derivative of lifted state with respect to x

dlift = jacobian(Basis,x);
jacfunName = [ 'jacobianLift_' , params.basisID , num2string(params.maxDegree) ];
matlabFunction(dlift, 'File', jacfunName, 'Vars', {zeta});

%% output variables  
params.Basis = Basis;    % symbolic vector of basis functions, p(x)
params.jacobianBasis = dlift;   % symbolic jacobian of the vector of basis monomials
params.N = N;   % dimension of basis (including the state itself)
params.Np = N + p;  % dimension of the lifted state
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xd = xd; % symbolic state delays
params.ud = ud; % symbolic input delays
params.liftHandle = liftFunName;    % name of the lifting function

end