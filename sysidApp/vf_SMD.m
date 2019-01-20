function [dx,y] = vf_SMD(t,x,u,m,b,k,params)
% DCmotor_vf: Ordinary differential equations for double pendulum.
%
%   author:  Alexander Erlich (alexander.erlich@gmail.com)
%
%   parameters:
%
%   t       Column vector of time points 
%   xdot    Solution array. Each row in xdot corresponds to the solution at a
%           time returned in the corresponding row of t.
%
%
%   ---------------------------------------------------------------------

dx = [0, 1; -k/m, -b/m] * x + [0, 1]' * u;

y = x;

end