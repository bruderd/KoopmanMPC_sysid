% main_MPC
%
%
%
%
%

%% load in system model
[model_file , model_path] = uigetfile('./models/');
model = load( [model_path , model_file] );

%% set MPC parameters
mpc = struct;

mpc.Ts      = model.params.Ts;      % sampling time must be the same as discrete model
mpc.tf      = 20;                     % total length of MPC simulation
mpc.tspan   = 0 : mpc.Ts : mpc.tf;  % time vector for MPC simulation
mpc.Np      = 100;                 % prediction horizon
mpc.x0      = zeros(model.params.n , 1);    % initial condition
mpc.nc      = 1;                % number of constraints


%% define reference trajectory

for i = 1 : length(mpc.tspan)
    ref = [ 0.08 * sin( mpc.tspan(i) ) ];   % sine wave
%     ref = [ 0.1 ];      % step input
    Yr( (i-1)*model.params.ny+1 : i*model.params.ny , 1 ) = ref;
end


%% define cost function matrices
% Cost function is defined: U'HU + ( z0'G + Yr'D )U

% A
N = size(model.A,1);
A = sparse( N*(mpc.Np+1) , N );
for i = 0 : mpc.Np
    A( (N*i + 1) : N*(i+1) , : ) = model.A^i ;
end

% B
Bheight = N*(mpc.Np+1);
Bcol = sparse( Bheight , size(model.B,2) );    % first column of B matrix
for i = 1 : mpc.Np
    Bcol( (N*i + 1) : N*(i+1) , : ) = model.A^(i-1) * model.B ;
end

Lshift = spdiags( ones( N*mpc.Np , 1 ) , -N , N*(mpc.Np+1) , N*(mpc.Np+1) );    % lower shift operator

Bwidth = size(model.B,2)*(mpc.Np);
B = spalloc( Bheight , Bwidth , floor(Bheight * Bwidth / 2) ); % initialze sparse B matrix
B(:,1) = Bcol;
for i = 2 : Bwidth
    B(:,i) = Lshift * B(:,i-1);
end

% C
C = kron( speye(mpc.Np+1) , model.C);

% Q
Q = kron( speye(mpc.Np+1) , eye(model.params.ny) * 1); % error magnitude penalty

% R
R = kron( speye(mpc.Np) , eye(model.params.p) * 0.0001);  % input magnitude penalty

% H, G, D
H = B' * C' * Q * C * B + R;
G = 2 * A' * C' * Q * C * B;
D = -2 * Q * C * B;

%% define constraint matrices 

nc = mpc.nc;

% F
F = sparse( nc * (mpc.Np+1) , size(B,2) );  % no constraints, all zeros

% E
E = sparse( nc * (mpc.Np+1) , size(B,1) );  % no constraints, all zeros

% L , M
L = F + E*B;
M = E*A;

%% save matrices in struct
mpc.H = H; mpc.G = G; mpc.D = D; mpc.L = L; mpc.M = M;


%% run MPC simulation

% identify system dynamics function
[vf_file , vf_path] = uigetfile('./fakeSystems/simDynamics');
fileparsed = split( vf_file , '.' );
system_dynamics = str2func( fileparsed{1} );

addpath(vf_path);
[t, x] = ode45( @(t,x) system_dynamics(x, get_MPCinput(t,x,Yr,mpc,model)), mpc.tspan(1 : end-mpc.Np), mpc.x0 );
rmpath(vf_path);
   
clear get_MPCinput; % clears the persistent variables from the get_MPCinput function

%% run open loop simulation

% % identify system dynamics function
% [vf_file , vf_path] = uigetfile('./fakeSystems/simDynamics');
% fileparsed = split( vf_file , '.' );
% system_dynamics = str2func( fileparsed{1} );
% 
% % solve for input over entire time horizon based on model
% u = gen_OLinput( 0 , mpc.x0 , Yr , mpc , model);
% 
% % simulate system
% addpath(vf_path);
% [t, x] = ode45( @(t,x) system_dynamics(x, get_OLinput(t,u,mpc,model)), mpc.tspan(1 : mpc.Np), mpc.x0 );
% rmpath(vf_path);


