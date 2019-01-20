% main_sysidApp_iter
% 
% Identify a nonlinear state space representation from data using the sysid
% App rather than the Koopman sysid method

% load the file of the Koopman model
model = load('C:\Users\danie\Documents\MATLAB\UM-Matlab\Koopman\ver8-Koopman_MPC\models\waves_192val_larm_sc09_191000pts_1del_Ts1_poly4_3.mat');

%% load data file
load('C:\Users\danie\Documents\MATLAB\UM-Matlab\Koopman\ver8-Koopman_MPC\dataFiles\larm_192val_16sid_sc09_191000pts_1del_Ts1.mat');
% data = load('C:\Users\danie\Documents\MATLAB\UM-Matlab\Koopman\ver8-Koopman_MPC\fakeSystems\simData\SMD_100s.mat');

%% construct iddata object
cd('..');
[zsysid_merged, zval_merged, zsysid, zval] = prep_iddata(data);
cd('sysidApp');
% zsysid_merged = iddata( data.x , data.u , 0.01 );
% zsysid_merged = iddata( data.alltrials.x , data.alltrials.u , model.params.Ts );


%% construct idnlgray object
Order = [ model.params.ny , model.params.p , 2 * model.params.n ];      % [Ny Nu Nx]
Parameters = {1e-6 * ones(model.params.N,1) , 1e-6 * ones(model.params.N,1)};   % initial parameter estimates
% Parameters = [1 1 1];
InitialStates = zeros( 2 * model.params.n , 1 ); % initial states
Ts = model.params.Ts;
% Ts = 0.01;
m = idnlgrey( 'vf_poly' , Order , Parameters , InitialStates , Ts , 'Name' , 'laser_polyModel');
% m = idnlgrey( 'vf_SMD' , Order , Parameters , [] , Ts , 'Name' , 'SMD');

%% merge 10 of the validation trials (20s of data to learn from total)

% initialize merged dataset
ztry = zval.z1;
for i = 2 : 10
   expID = ['z', num2str(i)];
   
   % merge all of the data sets into single multiexperiment object
   ztry = merge( ztry, zval.(expID) );
end

%% learn model
opt = nlgreyestOptions;
opt.Display = 'on';
% opt.SearchMethod  = 'fmincon';
nlmodel = nlgreyest( ztry , m , opt );

%% repeat until all of the validation data is used for training
numVals = numel(fieldnames(zval));
for i = 2 : floor( numvals / 10 )   % take them 10 at a time

    % the next 10 validation trials
    expID = ['z' , num2str(10*i+1)];
    ztry = zval.(expID);
    for j = 2:10
        expID = ['z' , num2str(10*i+j)];
        ztry = merge( ztry , zval.(expID) );
    end
    
    % intialize the new model from the old model
    Parameters = nlmodel.Parameters;
    InitialStates = nlmodel.InitialStates;
    
    % find new model
    nlmodel = nlgreyest( ztry , nlmodel , opt );
    
end
    

%% Check model accuracy

% estimate initial condition, don't just set to zero
compopt = compareOptions('InitialCondition','e'); 
for i = 1:4
m2.InitialStates(i).Fixed = false;
end

[y,fit,x0] = compare( zval_merged , nlmodel , compopt );
