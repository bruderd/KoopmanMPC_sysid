% gen_OLinput: function to determine open loop input over entire time horizon
function input = gen_OLinput( tnow , xnow , ref , mpc , model)

know = floor( tnow / mpc.Ts );    % current timestep

nd = model.params.nd;
if nd == 0
    zeta = xnow;
else
    zeta = [ xnow ; kron( ones(nd,1) , xnow ) ; zeros(nd * model.params.p , 1) ];   % assume state doesn't change during all time before 0
end

% lift the state
cd( 'liftingFunctions' );
lift = str2func( [ 'lift_' , model.params.systemName ] );
cd('..');
z = lift(zeta);   % lift x

Yr = ref( (know)*model.params.ny+1 : ((know+1) + mpc.Np)*model.params.ny , 1 ); % reference trajectory over prediction horizon
H = 2 * mpc.H;
f = ( z' * mpc.G + Yr' * mpc.D )';
A = mpc.L;
b = - mpc.M * z;
U = quadprog( H , f , A , b );     % use H, G etc here...

input = reshape(U , [model.params.p , mpc.Np]);     % inputs over whole time horizon [u1, u2, ...]

end