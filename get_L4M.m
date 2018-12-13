function L = get_L4M( A , B , koopData , params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Px = koopData.Px;
U = koopData.u;

L = zeros( size(Px,1) , params.N );
for i = 1 : size( Px , 1)    
    L(i,:) = ( A * Px(i,:)' + B * U(i,:)' )' ;        % with input
%     L(i,:) = ( A * Px(i,:)' )' ;        % without input
end

end

