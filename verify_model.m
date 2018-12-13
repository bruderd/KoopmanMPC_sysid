function verData = verify_model( tobs, xobs , uobs )
%verify_model: Simulates a model under real observed inputs and compares
%the simulation to what the system actually did.
%   Note: Can only use this with data that is properly sampled and scaled

%% initializations
% Load model and other system parameters into workspace
[modelFile, modelPath] = uigetfile('models');    % load model parameters
model = load([modelPath , filesep , modelFile]);    % load model file

Ts = model.params.Ts;    % sampling time of model

% initialize mpc data struct
verData = struct;
verData.T = [];
verData.U = [];
verData.Y = [];
verData.K = [];

params = model.params;    % model parameters

% get handle for the lifting function
modelName = modelFile(1:end-4); % remove the .mat file extension
cd( 'liftingFunctions' );
stateLift = str2func( [ 'lift_' , modelName ] );
cd( '..' );

% scale things down appropriately
xobs_sc = xobs ;% * diag( model.params.xScaleFactor );
uobs_sc = uobs ;% * diag( model.params.uScaleFactor );


%% real system

index0 = params.nd + 1;  % index of the first state
tspan = tobs(index0 : end) - tobs(index0);    % start simulation late so delays can be taken into account
[treal, xreal] = deal(tspan, xobs(index0 : end , :));


%% simulate the behavior of the learned lifted system

% set initial condition
x0 = xobs_sc(index0 , :)';
xd0 = reshape( flipud( xobs_sc(1 : params.nd , :) )' , [params.n * params.nd , 1] );
ud0 = reshape( flipud( uobs_sc(1 : params.nd , :) )' , [params.p * params.nd , 1] );
zeta0 = [x0; xd0; ud0];

%% simulate the behavior of the discrete linear system
xdis = zeros(length(tspan) , params.n);
xdis(1,:) = xobs_sc(index0 , :);
for i = 1 : length(tspan)-1
    if i == 1
        zetak = zeta0;
    else
        xk = xdis(i , :)';
        xdk = reshape( flipud( xdis( (i-params.nd) : i-1 , : ) )' , [params.n * params.nd , 1] );
        udk = reshape( flipud( uobs_sc(i-params.nd : i-1 , :) )' , [params.p * params.nd , 1] );
        zetak = [xk; xdk; udk];
    end
    psik = stateLift(zetak);
    psikp1 = model.Asim * psik + model.Bsim * uobs(i,:)';
    xdis(i+1,:) = ( model.C * psikp1 )';
end

%% quantify the error between real behavior and simulated behavior

%     % quantify error
%     terror = treal;
%     xerror = abs( xreal - xsysid );
% %     xerror = abs( xreal - xkoop );
%     xerrormax = max(max(xerror(:,1:ceil(params.n/2))));
%     % xerrormin = min(min(xerror(:,1:ceil(params.n/2))));
%     RMSE = sqrt( sum( (xreal - xsysid).^2 ) / length(terror) );
error = 0;

%% define outputs (and scale back up)

% simulated system
verData.T = tspan;
verData.Y = xdis * diag( model.params.xScaleFactor.^(-1) );
verData.U = uobs(index0 : end , :) * diag( model.params.uScaleFactor.^(-1) );

% real system
verData.Treal = tspan;
verData.Yreal = xreal * diag( model.params.xScaleFactor.^(-1) );
verData.Ureal = verData.U * diag( model.params.uScaleFactor.^(-1) );


%% plot the results

figure;

subplot(2,2,1)
hold on
p2 = plot( verData.Y(:,1) , verData.Y(:,2) );
p3 = plot( verData.Yreal(:,1) , verData.Yreal(:,2) );
p1 = plot( verData.Y(1,1) , verData.Y(1,2) , '*');
hold off
ylim([-8,8]);
xlim([-8,8]);
title('Trajectory of Real vs. Simulated')
legend([p2,p3] , {'Sim' , 'Real'} , 'Location' , 'northeastoutside')

subplot(2,2,2)
plot( verData.T , verData.U );
ylim([0,10]);
title( 'Control Inputs' )

subplot(2,2,3)
hold on
plot( verData.T , verData.Y(:,1) );
plot( verData.T , verData.Yreal(:,1) );
hold off
ylim([-8,8]);
ylabel('x-coordinate')
xlabel('time (s)')
legend({ 'Sim' , 'Real'} , 'Location' , 'northeastoutside')

subplot(2,2,4)
hold on
plot( verData.T , verData.Y(:,2) );
plot( verData.T , verData.Yreal(:,2) );
hold off
ylim([-8,8]);
ylabel('y-coordinate')
xlabel('time (s)')
legend({ 'Sim' , 'Real'} , 'Location' , 'northeastoutside')


end