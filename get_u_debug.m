function u = get_u_debug(t, x, valdata, valparams)
%get_u: Interpolates to estimate the value of the input at a specific t
% JUST USE FOR DEBUGGING, THIS FUNCTION IS NEVER ACTUALLY CALLED

u = interp1(valdata.t, valdata.u, t)';

end