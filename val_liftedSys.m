function [error, koopsim] = val_liftedSys( data, lifted )
%val_liftedSys: Validates the lifted system by comparing simulations to
%validation trials
%   Detailed explanation goes here

stateLift = str2func( lifted.params.liftHandle );   % identify handle of the lifting function

params = data.valparams;    % model parameters

error = struct;     % error results comparing real and koopman system
error.RMSE.total = 0; % initialize total error quantity
error.L1.total = 0;
error.dist.total = 0;
error.L1zero.total = 0;
error.L2zero.total = 0;
error.distzero.total = 0;

koopsim = struct;   % simulation results for koopman system

for j = 1 : params.numVals
    %% real system
     
    % isolate the jth validation trial
    valID = ['val', num2str(j)];
    valdata = data.(valID);
    
    index0 = params.nd + 1;  % index of the first state
    tspan = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
    [treal, xreal] = deal(tspan, valdata.x(index0 : end , :));
    
    
    %% simulate the behavior of learned state space system
%     
%     % set initial condition
%     x0 = valdata.x(index0 , :)';
% 
%     % simulate state space model
%     [tsysid, xsysid] = ode45(@(t,x) vf_koopman(x, get_delays(t, x, valdata, params) , get_u(t, x, valdata, params)), tspan, x0);
%     clear tpast xpast upast;    % clear persistent variables
% %     [tsysid, xsysid] = deal(treal, xreal);  % JUST A PLACEHOLDER. USE THIS IF WANT TO AVOID INTEGRATION ERROR
% 

    %% simulate the behavior of the learned lifted system
    
    % set initial condition
    x0 = valdata.x(index0 , :)';
    xd0 = reshape( flipud( valdata.x(1 : params.nd , :) )' , [params.n * params.nd , 1] );
    ud0 = reshape( flipud( valdata.u(1 : params.nd , :) )' , [params.p * params.nd , 1] );
    zeta0 = [x0; xd0; ud0];
    
    % simulate lifted linear model using ode45
%     [xss,tlifted,xlifted] = lsim(lifted.sys, valdata.u(index0 : end , :) , valdata.t(index0 : end , :) , stateLift(zeta0));
    
    % simulate the behavior of the discrete linear system (flow initial state forward, which is how MPC does it)
    xdis = zeros(length(tspan) , params.n);
    xdis(1,:) = valdata.x(index0 , :);
    psi0 = stateLift(zeta0);
    for i = 1 : length(tspan)-1
        if i == 1
            psik = psi0;
        else
            psik = psikp1;
        end
        psikp1 = lifted.A * psik + lifted.B * valdata.u(i,:)';
        xdis(i+1,:) = ( lifted.C * psikp1 )';
    end

%     % simulate the behavior of the discrete linear system (lift and drop at each step)
%     xdis = zeros(length(tspan) , params.n);
%     xdis(1,:) = valdata.x(index0 , :);
%     for i = 1 : length(tspan)-1
%         if i == 1
%             zetak = zeta0;
%         else
%             xk = xdis(i , :)';
%             xdk = reshape( flipud( xdis( (i-params.nd) : i-1 , : ) )' , [params.n * params.nd , 1] );
%             udk = reshape( flipud( valdata.u(i-params.nd : i-1 , :) )' , [params.p * params.nd , 1] );
%             zetak = [xk; xdk; udk];
%         end
%         psik = stateLift(zetak);
%         psikp1 = lifted.Asim * psik + lifted.Bsim * valdata.u(i,:)';
%         xdis(i+1,:) = ( lifted.C * psikp1 )';
%     end
    
    %% quantify the error between real behavior and simulated behavior
    
    % quantify L2 error
    terror = treal;
    error.RMSE.(valID) = sqrt( sum( (xreal - xdis).^2 ) / length(terror) );     % error for this trial
    error.RMSE.total = error.RMSE.total + error.RMSE.(valID);   % keep track of total error
    
    % quantify L1 error
    terror = treal;
    error.L1.(valID) = sum( abs(xreal - xdis) ) / length(terror);     % error for this trial
    error.L1.total = error.L1.total + error.L1.(valID);   % keep track of total error

    % quantify average cartesian distance error
    terror = treal;
    error.dist.(valID) = sum( sqrt( sum( (xreal - xdis).^2 , 2) ) ) / length(terror);     % error for this trial
    error.dist.total = error.dist.total + error.dist.(valID);   % keep track of total error
    
    % compute L2 error of the zero solution (if state remains zero whole time). This will be used for normalization purposes 
    error.L2zero.(valID) = sqrt( sum( (xreal).^2 ) / length(terror) );
    error.L2zero.total = error.L2zero.total + error.L2zero.(valID);   % keep track of total
    
    % compute L1 error of the zero solution (if state remains zero whole time). This will be used for normalization purposes 
    error.L1zero.(valID) = sum( abs(xreal) ) / length(terror);
    error.L1zero.total = error.L1zero.total + error.L1zero.(valID);   % keep track of total

    % compute average cartesian distance error of the zero solution (if state remains zero whole time). This will be used for normalization purposes 
    error.distzero.(valID) = sum( sqrt( sum( xreal.^2 , 2) ) ) / length(terror);
    error.distzero.total = error.distzero.total + error.distzero.(valID);   % keep track of total
    
    %% define outputs
%     error.(valID).terror = terror;
%     error.(valID).xerror = xerror;
%     error.(valID).RMSE = RMSE;
    
    koopsim.(valID).t = tspan;
%     koopsim.(valID).x = xsysid;       % should make this include both later...
%     koopsim.(valID).x = xss;
    koopsim.(valID).x = xdis;
    koopsim.(valID).u = valdata.u(index0 : end , :);
    
    
    %% plot the results
    
%     if params.ploton
%         figure
%         subplot(3,1,1)
%         plot(treal, xreal(:,1:ceil(params.n)))
%         title('Real system')
%         subplot(3,1,2)
%         plot(tsysid, xsysid(:,1:ceil(params.n)))
%         title('Identified system')
%         subplot(3,1,3)
%         hold on
%         plot(terror, xerror(:,1:ceil(params.n)))
%         plot(terror, xerrormax * ones(size(terror)), '--')
%         title('Error')
%         hold off
%     end

end

end


function u = get_u(t, x, valdata, params)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(valdata.t, valdata.u, t)';

end

function delays = get_delays(t, x, valdata, params)
    
    persistent tpast xpast upast
    
    % store the past values of t,x,u at ode solver steps
    if isempty(tpast)
        tpast = valdata.t(1 : params.nd+1);
        xpast = valdata.x(1 : params.nd+1 , :);
        upast = valdata.u(1 : params.nd+1 , :);
    else
        tpast = [tpast ; t];
        xpast = [xpast ; x'];
        upast = [upast ; get_u(t,x,valdata,params)'];
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
    xq_nd = xq_sofar(end-params.nd : end , :);
    uq_nd = uq_sofar(end-params.nd : end , :);
    
    % vectorize
    xd = reshape( flipud( xq_nd(1 : params.nd , :) )' , [params.n * params.nd , 1] );
    ud = reshape( flipud( uq_nd(1 : params.nd , :) )' , [params.p * params.nd , 1] );
    
    % set the output
    delays = [xd; ud];

end


