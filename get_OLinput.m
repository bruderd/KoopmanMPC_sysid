% get_OLinput: selects the current input from a vector of future inputs
function unow = get_OLinput( tnow , input , mpc , model)

know = floor( tnow / mpc.Ts ) + 1;    % current timestep

unow = input( : , know );

end