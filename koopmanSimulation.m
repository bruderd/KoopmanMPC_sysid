function [error, xkoop] = koopmanSimulation( valdata, valparams, koopman )
%koopmanValidation: Summary of this function goes here
%   Detailed explanation goes here

%% simulate the behavior of learned system

% only want to simulate for about 30 seconds
% i_30s = floor(30 / valparams.Ts);  % index of roughly 30 seconds in
% tspan = valdata.t(1 : i_30s);
tspan = valdata.t;

x0sim = valdata.x(1,:)'; % same initial state as validation data initial state

% simulated forward using the transpose of Koopman operator
xselector = [zeros(valparams.n,1), eye(valparams.n), zeros(valparams.n, valparams.N - valparams.n - 1)]; % matrix to extract state from lifted state
xkoop = zeros(length(tspan), valparams.n);
xkoop(1,:) = x0sim';
for i = 2 : length(tspan)
    ti = tspan(i);
    xnext = xselector * koopman.U' * polyLift( xkoop(i-1,:)' , get_u(ti, 0, valdata, valparams) );
    xkoop(i,:) = xnext';
end

% real system
[treal, xreal] = deal(tspan, valdata.x);


%% quantify the error between real behavior and simulated behavior

terror = treal;
xerror = abs( xreal - xkoop );
xerrormax = max(max(xerror(:,1:ceil(valparams.n/2))));
% xerrormin = min(min(xerror(:,1:ceil(valparams.n/2))));
RMSE = sqrt( sum( (xreal - xkoop).^2 ) / length(terror) );


%% plot the results

if valparams.ploton
    figure
    subplot(3,1,1)
    plot(treal, xreal(:,1:ceil(valparams.n/2)))
    title('Real system')
   subplot(3,1,2)
    plot(tspan, xkoop(:,1:ceil(valparams.n/2)))
    title('Koopman Transpose Sim')
    subplot(3,1,3)
    hold on
    plot(terror, xerror(:,1:ceil(valparams.n/2)))
    plot(terror, xerrormax * ones(size(terror)), '--')
    title('Error')
    hold off
end


% % animate the results
% animate_doublePendulum(sol_real, sol_sysid, valparams);

%% Define outputs

error = struct;
error.terror = terror;
error.xerror = xerror;
error.RMSE = RMSE;


end


function u = get_u(t, x, valdata, valparams)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(valdata.t, valdata.u, t)';

end