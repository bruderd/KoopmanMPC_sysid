% function to determine input at each time step
function unext = get_MPCinput( tnow , xnow , ref , mpc , model)

persistent k x u

know = floor( tnow / mpc.Ts );    % current timestep

% initialize values of persistent variables
if isempty(k)
    k = know;
    x = kron( ones(model.params.nd+1 , 1) , mpc.x0 );   % repeat initial condition for time steps that occur before k = 0
    u = kron( ones(model.params.nd+1 , 1) , zeros(model.params.p , 1) );    % assume input is zero before k = 0
end

% update values of persistent variables
if know == k(end)       % still at same timestep?
    unext = u(: , end);     % u takes same value 
else
    nd = model.params.nd;
    if nd == 0
        zeta = xnow;
    else
        zeta = [ xnow ; x( : , end-nd+1 : end ) ; u( : , end-nd+1 : end )];
    end
    
    % lift the state
    cd( 'liftingFunctions' );
    lift = str2func( [ 'lift_' , model.params.systemName ] );
    cd('..');
    z = lift(zeta);   % lift zeta

    Yr = ref( (know)*model.params.ny+1 : ((know+1) + mpc.Np)*model.params.ny , 1 ); % reference trajectory over prediction horizon
    H = 2 * mpc.H;
    f = ( z' * mpc.G + Yr' * mpc.D )';
    A = mpc.L;
    b = - mpc.M * z;
    U = quadprog( H , f , A , b );     % use H, G etc here...
    unext = U( 1 : model.params.p );    % the actual output of this function
    
    % store values of the state and input at the "beginning" of a timestep
    k = [k , know];     % vector of all previous timestamps
    x = [x , xnow];     % vector of all previous states
    u = [u , unext];
end


end