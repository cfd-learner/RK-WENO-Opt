function Amat = GetFDOperator( N )
%% GETFDOPERATOR Generates the matrix corresponding to the conservative
%                finite-difference method that maps the interface fluxes
%                to the cell-centered derivatives

Amat = zeros(N,N+1);
for i=1:N
    Amat(i,i)   = -1.0;
    Amat(i,i+1) =  1.0;
end

%% Done
end

