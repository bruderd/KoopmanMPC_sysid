function [ U , koopData ] = get_KoopmanConstGen( snapshotPairs , params )
%get_KoopmanConstGen: Find the best possible koopman operator given
%snapshot pairs using constraint generation to deal with large data sets.
%   Detailed explanation goes here

disp('Finding Koopman operator approximation...');

stateLift = str2func( params.liftHandle );   % identify handle of the lifting function

%% Extract snapshot pairs

% [x,y,u] = deal(snapshotPairs.x, snapshotPairs.y, snapshotPairs.u);
[x,y,u] = deal(snapshotPairs.zeta_x, snapshotPairs.zeta_y, snapshotPairs.u);

%% Build matrices

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.Np);
Py = zeros(length(x), params.Np);
for i = 1:length(x)
    psix = stateLift( x(i,:)' )';
    psiy = stateLift( y(i,:)' )';
    Px(i,:) = [ psix , u(i,:) ];
    Py(i,:) = [ psiy , zeros(1,p) ];     % exclude u from Py (could also use same u as Px
end

K = size(Px,1);
Np = size(Px,2);

%% Store useful data that can be used outside this function
koopData.Px = Px( : , 1 : params.N );   % only want state, not input
koopData.Py = Py( : , 1 : params.N );
koopData.x = snapshotPairs.x;
koopData.u = u;
koopData.zeta_x = snapshotPairs.zeta_x;

%% Solve for inital Koopman Operator wish subset of data points

% Build Apx sparsely with 10% of snapshotPairs
% Ktithe = min( floor(K/2) , 1000 );   % roughly 50% of total snapshotPairs, at most 1000 
Ktithe = K; % use all of the points

% Call function that solves QP problem
Uvec = solve_KoopmanQP(Px, Py, params);

U = reshape(Uvec, [Np,Np]);

%% check how well U works at each point (optional)
% dif = Px * U - Py;
% dif_x = dif( : , 1 : params.n);
    

%% Check which points the solution holds for, and repeat process as necessary

% optimal = false;
% satConst = zeros(K,1);      % logical array that stores which constraints are satisfied
% while ~optimal
%     unsatisfied = [];       % stores the indices of all the unsatisfied constraints
%     count = 0;
%     
%     % check if constraints are satisfied
%     for i = 1:K
%         if satConst(i,1) == 0         %~any(find(satisfied == i))   % only check points that weren't part of the last optimization problem
%             satConst(i,1) = all( Px(i,:) * U - Py(i,:) <= params.epsilon');        %check if the constraints are satisfied for this point
%             if satConst(i,1) == 0
%                 unsatisfied = [unsatisfied, i];     % store index of the unsatisfactory point
%                 count = count + 1;
%             end
%         end
%         % ensures we add at most 1000 new points to our optimization problem
%         if count > 1000
%             break;
%         end
%     end
%     
%     % print progress and check if all (most) of the constraints are satisfied
%     progress = 100 * sum(satConst)/K;
%     disp(['Epsilon = ', num2str(mean(params.epsilon))]);     % print the average value of epsilon
%     disp([num2str(progress), '% of constraints satisfied']);    % print percentage of constraints satisfied
%     if (sum(satConst) > params.percSat*K)
%         optimal = true;
%         break;
%     end
%     
%     % Construct new A and b matrices
%     Kunsat = length(unsatisfied);   % number of points where constraints unsatisfied
%     row = zeros(Kunsat*N^2,1);
%     col = zeros(Kunsat*N^2,1);
%     val = zeros(Kunsat*N^2,1);
%     for i = 1 : Kunsat
%         for j = 1 : N
%             index = N*(i-1) + j;
%             row(1+(index-1)*N : index*N) = index;
%             col(1+(index-1)*N : index*N) = mod(index-1,N) * N + 1 : mod(index-1,N) * N + N;
%             val((1+(index-1)*N : index*N)) = Px( unsatisfied(i) , : );
%         end
%     end
%     Apx_add = sparse(row, col, val, Kunsat*N, N*N);
%     bpy_add = reshape( Py(unsatisfied,:)' , [Kunsat*N,1]);
%     Apx = [Apx; Apx_add];
%     bpy = [bpy; bpy_add];
%     
%     % Call function to solve QP
%     Uvec = solve_KoopmanQP(Px, Py, params);
%     U = reshape(Uvec, [N,N]); 
%     
% end

disp('Done');

end

