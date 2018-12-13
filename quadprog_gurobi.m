function [ x , results ] = quadprog_gurobi( H, f, A, b )
%quadprog: Solves a quadratic program using Gurobi
%   Detailed explanation goes here

%% solve QP using gurobi

model.Q = sparse(0.5 * H);
model.A = sparse(A);
model.obj = f';
model.rhs = b';
model.sense = '<';

% solve using gurobi
results = gurobi(model);

% extract the solution from the results struct
x = results.x;

end

