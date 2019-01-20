function [dx,y] = vf_poly(t,x,u,cf1,cf2,params)
%vf_poly: nonlinear (polynomial) state-space dynamics for the flaccy-laser system  
%
% The output variables are:
% 
%     dx — Represents the right side(s) of the state-space equation(s). A column vector with Nx entries. For static models, dx=[].
% 
%     For discrete-time models. dx is the value of the states at the next time step x(t+Ts).
% 
%     For continuous-time models. dx is the state derivatives at time t, or dxdt.
% 
%     y — Represents the right side(s) of the output equation(s). A column vector with Ny entries.
% 
% The file inputs are:
% 
%     t — Current time.
% 
%     x — State vector at time t. For static models, equals [].
% 
%     u — Input vector at time t. For time-series models, equals [].
% 
%     cf1,cf2, ...,cfN — Parameters, which can be real scalars, column vectors or two-dimensional matrices. N is the number of parameter objects. For scalar parameters, N is the total number of parameter elements.
% 
%     params — Contains auxiliary variables that might be required for updating the constants in the state equations. 

%% dynamics

cd(['..' , filesep , 'liftingFunctions']);
% psi = str2func(['lift_' , params.systemName]);    % the lifted state, as a function
psi = str2func('lift_waves_192val_larm_sc09_191000pts_1del_Ts1_poly4_3');    % the lifted state, as a function
cd(['..' , filesep , 'sysidApp']);

%% limit the number of coefficients needed
vec1 = [cf1(1:119)' , zeros(1,210) , cf1(120)]; % total must add up to 330
vec2 = [cf2(1:119)' , zeros(1,210) , cf2(120)]; % total must add up to 330

%% define dx

dx(1,1) = vec1 * psi([x ; u']);
dx(2,1) = vec2 * psi([x ; u']);
dx(3,1) = x(1);
dx(4,1) = x(2);


%% output
% y = [eye(params.n) , zeros(params.n,params.n)] * x;
y = [zeros(2,2) , eye(2)] * x;
% y = x;