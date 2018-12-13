function R = get_R4M( L , Cz , koopData , params)
%UNTITLED2 Summary of this function goes here
%   FYI: this function calls 'stateLift' so must be in the same folder as it

Py = koopData.Py;

R = zeros( size(L,1) , params.N );
for i = 1 : size( L , 1 )  
%     R(i,:) = ( stateLift( Cz * L(1,:)' ) )' ;
     R(i,:) = Py(i,:) ;  % debug. This should make M the identity
end

end

