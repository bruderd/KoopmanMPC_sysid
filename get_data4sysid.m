function [ data4sysid ] = get_data4sysid( data , koopsim, params )
%get_data4sysid: Converts trial, validation, and simulation data to iddata
%format so that it can be compared to outputs of systems from matlab sysid
%toolbox

% convert data to a format matlabs sysid toolbox can use
[zsysid_merged, zval_merged, zsysid, zval] = prep_iddata(data);

% save in struct for output
data4sysid = struct;
data4sysid.sysid_merged = zsysid_merged;
data4sysid.val_merged = zval_merged;
data4sysid.val = zval;
data4sysid.sysid = zsysid;
for k = 1:params.numVals
    valID = ['val', num2str(k)];
    zID = ['z', num2str(k)];
    data4sysid.valkoop.(zID) = iddata(koopsim.(valID).x, koopsim.(valID).u, data.valparams.Ts, 'Name', 'Koopman');
end

% show comparison of Koopman system verses ground truth
if params.ploton
    for k = 1: params.numVals
        zID = ['z', num2str(k)];
        figure
        compare(data4sysid.val.(zID), data4sysid.valkoop.(zID));
        
        % change y-axis of figure
        fig = gcf;
        allaxes = findall(fig, 'type', 'axes');
        allaxes(5).YLim = [-params.scale , params.scale];
        allaxes(6).YLim = [-params.scale , params.scale];
%         ax = gca;
%         ax.YLim = [-params.scale , params.scale];    
    end
end

end

