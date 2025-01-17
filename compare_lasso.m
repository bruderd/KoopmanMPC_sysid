function results = compare_lasso( basisID , maxDegree , lassoParams )
%compare_lasso
%   Learns a model from data for a bunch of different values of the lassoParams
%   Used to compare the performance and sparsity of the models generated
%
%   Performs linear system identification of nonlinear systems using a
%   lifting technique based on the Koopman operator projected onto a finite
%   monomial basis.
%
%   INPUTS:
%       basisID - string that spedifies type of basis ('fourier' or 'poly' or 'fourier_sparser' or 'thinplate' or 'gaussian')
%       maxDegree - maximum "degree" of the basis elements
%       lassoParams - a vector containing all the differenta values for the
%           lasso optimizaton parameters that should be tried.
%
%   OUTPUTS:
%       lifted - a struct containing the lifted linear model matrices
%       error - a struct containing the error info for all of the
%       validation trials

addpath('/home/bruderd/gurobi/gurobi810/linux64/matlab');

%% Define system parameters (USER EDIT SECTION)

if ~exist('params' ,'var')  % recycle struct from previous run 
    params = struct;
end
% params.getData = 'manSim2_1val_1sid_sc09_100000pts_1del_Ts02.mat';            % (name of the data file) , for simulated 2d manipulator / double pendulum
params.getData = 'larm_10val100s_16sid_sc09_191000pts_1del_Ts1.mat';            % (name of the data file) , for laser flaccy system
% params.getData = 'larm_OL50pc_100val10s_16sid_sc09_191000pts_1del_Ts1.mat';            % (name of the data file) , for laser flaccy system
params.basisID = basisID;   % ('fourier' or 'poly' or 'fourier_sparser' or 'thinplate' or 'gaussian')

% parameters for reading in data (these affect how shapshot pairs built from raw data).
params.numTrials        = 16;        % numer of sysid trials
params.numVals          = 10;        % number of validation trials
params.Ts               = 0.1;     % sampling period
params.K                = 191000;     % numer of snapshotPairs to take
params.numericalDerivs  = false;    % choose whether or not to take numerical derivatives of states (boolean)
params.scale            = 0.9;      % scale down all state to be in range [-scale , scale]
params.nd               = 1;        % number of delays to include in the snapshot pairs

% params.systemName          = ['waves_manSim2_1val_1sid_sc09_100000pts_1del_Ts02_' , basisID , num2str(maxDegree) ];  % name of current system , for simulated 2d manipulator / double pendulum
params.systemName          = ['waves_larm_10val100s_16sid_sc09_191000pts_1del_Ts1_' , basisID , num2str(maxDegree) ];  % name of current system , for laser flaccy system
% params.systemName          = ['waves_larm_OL50pc_100val10s_16sid_sc09_191000pts_1del_Ts1_' , basisID , num2str(maxDegree) ];  % name of current system , for laser flaccy system, with outliers
params.filterWindow        = floor( [1/params.Ts, 1/params.Ts] );  % if taking numerical derivatives, specifies the moving mean window before and after derivatives taken.

% Koopman Sysid parameters
params.n = 2;   % dimension of state space (including state derivatives)
params.p = 3;   % dimension of input
params.ny = 2;  % dimension of output
params.naug = params.n + params.p; % dimension of augmented state (DNE)
params.nzeta = params.n + params.nd * (params.naug);    % dimensinon of zeta (DNE)

% select maximum "degree" for basis elements
params.maxDegree = maxDegree;   % maximum degree of vector field monomial basis

% only do this if the Basis is not already defined. Will need to clear before running with a different basis or maxDegree
if ~isfield(params , 'Basis')   
    % define lifting function and basis
    disp('Defining basis of observables...')
    if strcmp(params.basisID, 'fourier')
        params = def_fourierLift(params);  % creates fourier lifting function, fourierLift;
    elseif strcmp(params.basisID, 'poly')
        params = def_polyLift(params);  % creates polynomial lifting function, polyLift
    elseif strcmp(params.basisID, 'fourier_sparser')
        params = def_fourierLift_sparser(params);
    elseif strcmp(params.basisID, 'thinplate')
        params = def_thinplateLift(params);
    elseif strcmp(params.basisID, 'gaussian')
        params = def_gaussianLift(params);
    end
    disp('Done.')
end


% Koopman sysid tuning parameters
params.t        = 10 * params.N; % penalty on model complexity
params.epsilon  = 1; % model accuracy tolerance (larger value = less accurate)
params.percSat  = 0.75;  % percentage of snapshot pairs that must satisfy accuracy tolerance

% output parameters
params.validateon          = true;  % boolean to decide whether to validate the model
params.ploton              = false;  % boolean to turn validation/error/compare plot on or off
params.compareon           = false;  % boolean to decide whether to convert to iddata and compare to validation data with Matlab compare function


%% Get training data

% Load in the training data file
disp('Retreiving training data...');
data_path = 'dataFiles';
matcontents = load([ data_path , filesep ,  params.getData ]); % must be a .mat file
data = matcontents.data;
for fn = fieldnames(params)'    % update params but remember the ScaleFactor fields
    data.valparams.(fn{1}) = params.(fn{1});
end
params = data.valparams;
disp('Done.')

% Get the snapshot pairs
disp([ 'Extracting random sampling of ', num2str(params.K), ' snapshot pairs from data...'  ])
some_snapshotPairs = get_randsnapshotPairs(params.K, data.snapshotPairs);
disp('Done.')

%% Learn the approximate Koopman operator and corresponding NL system
results = struct;
results.lassoParams = lassoParams;

for i = 1 : length(lassoParams)
    disp(['Learning model number ' , num2str(i) , ' of ' , num2str(length(lassoParams)) , '.' ] );
    
    params.t = lassoParams(i) * params.N;
    [model , current_error] = learn_koopmanModel( data, some_snapshotPairs , params );   % learn new model
    
    % save relevant parts of the model in output struct
    if length(lassoParams) <= 100   % only save the full model each time if problem is small enough
        results.model(i) = model;
        results.error(i) = current_error;
    end
    results.density.A(i) = nnz( abs( model.Asim ) > 1e-6 ) / numel(model.Asim);  % element considered 0 if magnitude less than or equal to 1e-6
    results.density.MA(i) = nnz( abs( model.A ) > 1e-6 ) / numel(model.A);
    results.totnerror.L2(i,:) = current_error.RMSE.total ./ current_error.L2zero.total; % note these are 1x2
    results.totnerror.L1(i,:) = current_error.L1.total ./ current_error.L1zero.total;
    results.totnerror.dist(i,:) = current_error.dist.total ./ current_error.distzero.total; % note this is 1x1
    
%     % keep this model only if the error is lower than the last one
%     if norm( current_error.RMSE.total ) < norm( error.RMSE.total )  % want average of x and y error to be better
%         lifted = model;
%         error = current_error;
%         lifted.error = error;   % save the error information as part of the model
%     end
end

results.params = params;    % save the params so I can check them later

%% Convert to iddata and compare to validation data
% if params.compareon
%     disp('Converting to iddata format...');
%     data4sysid = get_data4sysid( data , koopsim , params );
%     disp('Done.')
% end

%% Save the results file

% save datafile without overwriting previous files with same name
[unique_fname, change_detect] = auto_rename(['lasso', filesep, params.systemName, '.mat'], '_0');
save(unique_fname, '-struct', 'results');
