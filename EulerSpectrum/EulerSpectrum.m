%% Script to compute the pseudo-spectrum of the finite-difference
%  discretization operators.

clear all;
close all;
addpath('./eigtoollib');


%% Set up some parameters
N = 200;
method = 'crweno5';
boundary = 'aperiodic';
mach = 1.0/3.0;

%% Compute the discretization Jacobian
fprintf('Computing discretization Jacobian.\n');
A = GetInterpOperator(N,method,boundary);
B = GetFDOperator(N);
C = -B*A;

%% Compute spectrum
fprintf('Computing spectrum.\n');
lambda = eig(C);
lambda_1 = abs(mach)   * lambda; % |u|
lambda_2 = abs(1-mach) * lambda; % |u-a|
lambda_3 = abs(1+mach) * lambda; % |u+a|
lambda = [lambda_1;lambda_2;lambda_3];

%% Plot
figure(1);
plot(real(lambda),imag(lambda),'bo');
axis equal;
grid on;
hold on;

%% Save to file
filename = ['eigenvalues_',method,'.dat'];
fid = fopen(filename,'w');
for i = 1:size(lambda,1)
    fprintf(fid,'%+1.16e %+1.16e\n',real(lambda(i)),imag(lambda(i)));
end
fclose(fid);