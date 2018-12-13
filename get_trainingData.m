function [ data, some_snapshotPairs , params ] = get_trainingData( params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Load from data file or read data from raw measurements
disp('Retreiving training data...');

getData = params.getData;

if strcmp(getData, 'exp')
    addpath('generateData');
    [data , params] = gen_data_fromExp( params );
    rmpath('generateData')
elseif strcmp(getData, 'file')
    % Prompt user to identify data file
    [data_file,data_path] = uigetfile;
    matcontents = load([data_path, data_file]); % must be a .mat file
    data = matcontents.data;
%     data.valparams = params;    % makes sure that the params match the current problem
    for fn = fieldnames(params)'    % update params but remember the ScaleFactor fields
        data.valparams.(fn{1}) = params.(fn{1});
    end
    params = data.valparams;
end

disp('Done.')

%% take random subset of snapshot pairs
disp([ 'Extracting random sampling of ', num2str(params.K), ' snapshot pairs from data...'  ])

some_snapshotPairs = get_randsnapshotPairs(params.K, data.snapshotPairs);

disp('Done.')

end

