% main_sysidApp
% 
% Identify a nonlinear state space representation from data using the sysid
% App rather than the Koopman sysid method

% load the file of the Koopman model
cd([ '..' , filesep , 'models' ]);
model = load('waves_192val_larm_sc09_191000pts_1del_Ts1_poly4_3.mat');
cd([ '..' , filesep , 'sysidApp' ]);

%% load data file
cd([ '..' , filesep , 'dataFiles' ]);
load('larm_192val_16sid_sc09_191000pts_1del_Ts1.mat');
cd([ '..' , filesep , 'sysidApp' ]);


%% construct iddata object
% cd('..');
[zsysid_merged, zval_merged, zsysid, zval] = prep_iddata_allscaled(data);
% cd('sysidApp');
% zsysid_merged = iddata( data.x , data.u , 0.01 );
% zsysid_merged = iddata( data.alltrials.x , data.alltrials.u , model.params.Ts );

%% construct idnlgray object
Order = [ model.params.ny , model.params.p , 2 * model.params.n ];      % [Ny Nu Nx]
Parameters = {1e-6 * ones(36,1) , 1e-6 * ones(36,1)};   % initial parameter estimates
% Parameters = {1e-6 * ones(120,1) , 1e-6 * ones(120,1)};   % initial parameter estimates
% Parameters = {1e-6 * ones(model.params.N,1) , 1e-6 * ones(model.params.N,1)};   % initial parameter estimates
InitialStates = zeros( 2 * model.params.n , 1 ); % initial states
Ts = model.params.Ts;
m = idnlgrey( 'vf_poly' , Order , Parameters , InitialStates , Ts , 'Name' , 'laser_polyModel');


%% merge 100 of the validation trials (200s of data to learn from total)

% initialize merged dataset
ztry = zval.z1(1:19);   % just the first 19 samples
for i = 2 : 100
   expID = ['z', num2str(i)];
   
   % merge all of the data sets into single multiexperiment object
   ztry = merge( ztry, zval.(expID)(1:19) );    % just the first 19 samples
end

%% learn model
opt = nlgreyestOptions;
opt.Display = 'on';
% opt.SearchMethod  = 'fmincon';
nlmodel = nlgreyest( ztry , m , opt );

%% Check model accuracy
% 
% % estimate initial condition, don't just set to zero
% compopt = compareOptions('InitialCondition','e'); 
% for i = 1:4
% nlmodel.InitialStates(i).Fixed = false;
% end
% [y,fit,x0] = compare( zval_merged , nlmodel , compopt );

%% Save NLGREY model
save( ['nlModels' , filesep , 'laserModel_poly2_allsyids.mat'] , 'nlmodel' );