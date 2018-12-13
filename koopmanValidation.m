function [error, koopsim] = koopmanValidation( data, valparams, statespace, lifted )
%koopmanValidation: Summary of this function goes here
%   Detailed explanation goes here

error = struct;     % error results comparing real and koopman system
koopsim = struct;   % simulation results for koopman system

for j = 1 : valparams.numVals
    %% real system
     
    % isolate the jth validation trial
    valID = ['val', num2str(j)];
    valdata = data.(valID);
    
    index0 = valparams.nd + 1;  % index of the first state
    tspan = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
    [treal, xreal] = deal(tspan, valdata.x(index0 : end , :));
    
    
    %% simulate the behavior of learned state space system
%     
%     % set initial condition
%     x0 = valdata.x(index0 , :)';
% 
%     % simulate state space model
%     [tsysid, xsysid] = ode45(@(t,x) vf_koopman(x, get_delays(t, x, valdata, valparams) , get_u(t, x, valdata, valparams)), tspan, x0);
%     clear tpast xpast upast;    % clear persistent variables
% %     [tsysid, xsysid] = deal(treal, xreal);  % JUST A PLACEHOLDER. USE THIS IF WANT TO AVOID INTEGRATION ERROR
% 

    %% simulate the behavior of the learned lifted system
    
    % set initial condition
    x0 = valdata.x(index0 , :)';
    xd0 = reshape( flipud( valdata.x(1 : valparams.nd , :) )' , [valparams.n * valparams.nd , 1] );
    ud0 = reshape( flipud( valdata.u(1 : valparams.nd , :) )' , [valparams.p * valparams.nd , 1] );
    zeta0 = [x0; xd0; ud0];
    
    % simulate lifted linear model
    [xss,tlifted,xlifted] = lsim(lifted.sys, valdata.u(index0 : end , :) , valdata.t(index0 : end , :) , stateLift(zeta0));
    
    
    %% quantify the error between real behavior and simulated behavior
    
%     % quantify error
%     terror = treal;
%     xerror = abs( xreal - xsysid );
% %     xerror = abs( xreal - xkoop );
%     xerrormax = max(max(xerror(:,1:ceil(valparams.n/2))));
%     % xerrormin = min(min(xerror(:,1:ceil(valparams.n/2))));
%     RMSE = sqrt( sum( (xreal - xsysid).^2 ) / length(terror) );
    error = 0;
    
    %% define outputs
%     error.(valID).terror = terror;
%     error.(valID).xerror = xerror;
%     error.(valID).RMSE = RMSE;
    
    koopsim.(valID).t = tspan;
%     koopsim.(valID).x = xsysid;       % should make this include both later...
    koopsim.(valID).x = xss;
    koopsim.(valID).u = valdata.u(index0 : end , :);
    
    
    %% plot the results
    
%     if valparams.ploton
%         figure
%         subplot(3,1,1)
%         plot(treal, xreal(:,1:ceil(valparams.n)))
%         title('Real system')
%         subplot(3,1,2)
%         plot(tsysid, xsysid(:,1:ceil(valparams.n)))
%         title('Identified system')
%         subplot(3,1,3)
%         hold on
%         plot(terror, xerror(:,1:ceil(valparams.n)))
%         plot(terror, xerrormax * ones(size(terror)), '--')
%         title('Error')
%         hold off
%     end

end

end


function u = get_u(t, x, valdata, valparams)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(valdata.t, valdata.u, t)';

end

function delays = get_delays(t, x, valdata, valparams)
    
    persistent tpast xpast upast
    
    % store the past values of t,x,u at ode solver steps
    if isempty(tpast)
        tpast = valdata.t(1 : valparams.nd+1);
        xpast = valdata.x(1 : valparams.nd+1 , :);
        upast = valdata.u(1 : valparams.nd+1 , :);
    else
        tpast = [tpast ; t];
        xpast = [xpast ; x'];
        upast = [upast ; get_u(t,x,valdata,valparams)'];
    end
    
    % remove repeated values caused by solver taking small steps
    [tpast_unq, ind_unq] = unique(tpast);
    xpast_unq = xpast(ind_unq, :);
    upast_unq = upast(ind_unq, :);
    
    % interpolate to find values at valid sampling points
    xq = interp1(tpast_unq, xpast_unq, valdata.t, 'linear', 0); % return 0 outside range
    uq = interp1(tpast_unq, upast_unq, valdata.t, 'linear', 0); % return 0 outside range
    
    % exclude values outsite the interpolation range
    xq_sofar = xq( ( xq(:,1) ~= 0 ) , :);
    uq_sofar = uq( ( xq(:,1) ~= 0 ) , :);
    
    % just take the last nd values
    xq_nd = xq_sofar(end-valparams.nd : end , :);
    uq_nd = uq_sofar(end-valparams.nd : end , :);
    
    % vectorize
    xd = reshape( flipud( xq_nd(1 : valparams.nd , :) )' , [valparams.n * valparams.nd , 1] );
    ud = reshape( flipud( uq_nd(1 : valparams.nd , :) )' , [valparams.p * valparams.nd , 1] );
    
    % set the output
    delays = [xd; ud];

end


