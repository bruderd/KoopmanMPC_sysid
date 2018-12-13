function data = gen_data_fromExp_xonly( params )
%genData_fromExp: Generate system data and snapshot pairs from experimental
%data
%   Detailed explanation goes here

% initialize output struct
data = struct;
alltrials = struct;

%% Read in data from all trials

num = params.numTrials;
trialCount = 0;        % trial counter

alltrials.t = []; alltrials.y = []; alltrials.u = []; alltrials.x = [];
x = []; y = [];
for i = 1:num
    trialCount = trialCount + 1;    % increment trial counter
    
    % generate data from one simulation
    disp(['Please select .mat file corresponding to trial number ', num2str(i), '...']);
    trialData = get_data_xonly(params);
    
    % append this data to the "alltrials" field of data
    alltrials.t = [alltrials.t; trialData.t];     % time vector
    alltrials.y = [alltrials.y; trialData.y];     % state "measurements"
    alltrials.u = [alltrials.u; trialData.u];     % input
    alltrials.x = [alltrials.x; trialData.x];     % actual state
    
    xk = zeros(length(trialData.t), size(trialData.y,2) + size(trialData.u,2));
    yk = zeros(length(trialData.t), size(trialData.y,2) + size(trialData.u,2));
    for j = 1:length(trialData.t)-1
        xk(j,:) = [ trialData.y(j,:), trialData.u(j,:) ];
        yk(j,:) = [ trialData.y(j+1,:), trialData.u(j,:) ];
    end
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
    
    % save this trial data to the output struct
    trialID = ['trial', num2str(trialCount)];
    data.(trialID) = trialData;
    
end

% define snapshotPairs struct
snapshotPairs = struct;
snapshotPairs.x = x;
snapshotPairs.y = y;

%% Read in validation data set

disp(['Please select .mat file corresponding to validation trial...']);
validation = get_data_xonly(params);


%% Define output

data.alltrials = alltrials;     % saves data from all trials as a single timeseries
data.snapshotPairs = snapshotPairs;
data.validation = validation;   % trial that can be used for model validation
data.valparams = params;   % saves params used for validation so we can remember

%% save datafile without overwriting previous files with same name
% SaveWithNumber(['dataFiles', filesep, params.systemName, '.mat'], data);
[unique_fname, change_detect] = auto_rename(['dataFiles', filesep, params.systemName, '.mat'], '0');
save(unique_fname, 'data');

end