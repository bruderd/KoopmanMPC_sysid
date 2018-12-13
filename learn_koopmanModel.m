function [ model , error ] = learn_koopmanModel( snapshotPairs , params )
%learn_koopmanModel: learns a koopman model and quantifies the error of
%that model over all validation trials
%   Detailed explanation goes here

[U , koopData ] = get_KoopmanConstGen( snapshotPairs, params );
model          = sysid_liftedSys( U , params , koopData );

%% Simulate the results and compare to validation trial(s)
if params.validateon
    disp('Comparing to validation data set...');
    [error, koopsim] = val_liftedSys(data, model);
    disp('Done.')
end


end

