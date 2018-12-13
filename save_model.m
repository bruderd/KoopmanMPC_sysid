function [unique_fname] = save_model( system )
%save_model: saves a learned model as a .mat file
%   Detailed explanation goes here

systemName = system.params.systemName;

%% save datafile without overwriting previous files with same name
[unique_fname, change_detect] = auto_rename(['models', filesep, systemName, '.mat'], '-0');
save(unique_fname, '-struct', 'system');

end

