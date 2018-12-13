function some_snapshotPairs = get_randsnapshotPairs( num , snapshotPairs )
%get_randsnapshotPairs: Randomly extracts num snapshot pairs
%   Detailed explanation goes here

some_snapshotPairs = struct;

totalPairs = length(snapshotPairs.x);

s = RandStream('mlfg6331_64'); 
index = datasample(s , 1:totalPairs , num , 'Replace' , false);
% index = randi(totalPairs, num, 1);

some_snapshotPairs.x = snapshotPairs.x(index,:);
some_snapshotPairs.y = snapshotPairs.y(index,:);
some_snapshotPairs.u = snapshotPairs.u(index,:);
some_snapshotPairs.zeta_x = snapshotPairs.zeta_x(index,:);
some_snapshotPairs.zeta_y = snapshotPairs.zeta_y(index,:);

end

