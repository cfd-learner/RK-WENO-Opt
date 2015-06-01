function Amat = GetInterpOperator( N,method,boundary )
%% GETINTERPOPERATOR Returns the matrix corresponding to an 
%                    interpolation method that computes the 
%                    interface fluxes from cell-centered fluxes

if (strcmp(strtrim(method),'1'))
    Amat = FirstOrderUpwind(N,boundary);
elseif (strcmp(strtrim(method),'weno5'))
    Amat = WENO5(N,boundary);
elseif (strcmp(strtrim(method),'crweno5'))
    Amat = CRWENO5(N,boundary);
else
    Amat = zeros(N+1,N);
end


%% Done
end

function Amat = FirstOrderUpwind( N,boundary )
%% GETINTERPOPERATOR Returns the matrix corresponding to the WENO5 scheme

Amat = zeros(N+1,N);

for i = 2:N+1
    Amat(i,i-1) = 1.0;
end

if (strcmp(strtrim(boundary),'periodic'))
    Amat(1,N) = 1.0;
end

%% Done
end

function Amat = WENO5( N,boundary )
%% GETINTERPOPERATOR Returns the matrix corresponding to the WENO5 scheme

c = [1.0/30.0, -13.0/60.0, 47.0/60.0, 27.0/60.0, -1.0/20.0];

Amat = zeros(N+1,N);

Amat(1,1:2) = c(4:5);
Amat(2,1:3) = c(3:5);
Amat(3,1:4) = c(2:5);

for i = 4:N-1
    Amat(i,i-3:i+1) = c(:);
end

Amat(N,N-3:N) = c(1:4);
Amat(N+1,N-2:N) = c(1:3);

if (strcmp(strtrim(boundary),'periodic'))
    Amat(1,N-2:N) = c(1:3);
    Amat(2,N-1:N) = c(1:2);
    Amat(3,N) = c(1);
    Amat(N,1) = c(5);
    Amat(N+1,1:2) = c(4:5);
end

%% Done
end

function Amat = CRWENO5( N,boundary )
%% GETINTERPOPERATOR Returns the matrix corresponding to the CRWENO5 scheme

lhs = [0.3, 0.6, 0.1];
rhs = [1.0/30.0, 19.0/30.0, 1.0/3.0];

L = zeros(N+1,N+1);
R = zeros(N+1,N);

if (strcmp(strtrim(boundary),'periodic'))
    
    L(1,1:2) = lhs(2:3);
    L(1,N+1) = lhs(1);
    for i = 2:N
        L(i,i-1:i+1) = lhs(:);
    end
    L(N+1,N:N+1) = lhs(1:2);
    L(N+1,1) = lhs(3);
    
    R(1,1) = rhs(3);
    R(1,N-1:N) = rhs(1:2);
    R(2,1:2) = rhs(2:3);
    R(2,N) = rhs(1);
    for i=3:N
        R(i,i-2:i) = rhs(:);
    end
    R(N+1,N-1:N) = rhs(1:2);
    R(N+1,1) = rhs(3);
    
else
    
    L(1,1) = 1.0;
    for i = 2:N
        L(i,i-1:i+1) = lhs(:);
    end
    L(N+1,N+1) = 1.0;
    
    R(1,1:2) = [27.0/60.0, -1.0/20.0]; % WENO5
    R(2,1:2) = rhs(2:3);
    for i=3:N
        R(i,i-2:i) = rhs(:);
    end
    R(N+1,N-2:N) = [1.0/30.0, -13.0/60.0, 47.0/60.0]; %WENO5
    
end

Amat = inv(L)*R; %#ok<MINV>

%% Done
end