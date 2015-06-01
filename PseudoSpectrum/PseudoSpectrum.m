%% Script to compute the pseudo-spectrum of the finite-difference
%  discretization operators.

clear all;
close all;
addpath('./eigtoollib');


%% Set up some parameters
N = 200;
method = 'weno5';
boundary = 'aperiodic';

%% Compute the discretization Jacobian
fprintf('Computing discretization Jacobian.\n');
A = GetInterpOperator(N,method,boundary);
B = GetFDOperator(N);
C = -B*A;

%% Compute and plot the spectrum
fprintf('Computing spectrum.\n');
lambda = eig(C);
figure(1);
plot(real(lambda),imag(lambda),'bo');
axis equal;
grid on;
hold on;

%% Compute and plot the pseudo-spectrum
fprintf('Computing pseudo-spectrum.\n');
eps_range = [-16,-14,-12,-10,-8,-6,-4];
opts.levels = eps_range;
opts.npts = 100;
opts.no_graphics = 1;
opts.no_waitbar = 1;
[x,y,sigs] = eigtool(C,opts);
figure(1);
[Contour,h] = contour(x,y,log10(sigs),'LevelList',eps_range);
caxis([min(eps_range),max(eps_range)]);
colorbar;

%% Save the pseudo-spectrum
fprintf('Saving pseudo-spectrum.\n');
max = size(Contour,2);
i = 1;
filename = strcat('pseudospectrum_',strtrim(method),'.dat');
fid = fopen(filename,'w');
while true

    if (i > max)
        break;
    end

    Level = Contour(1,i);
    N = Contour(2,i);
    xcoord = zeros(N,1);
    ycoord = zeros(N,1);
    for j = 1:N
        xcoord(j) = Contour(1,i+j);
        ycoord(j) = Contour(2,i+j);
        fprintf(fid,'%+1.16e %+1.16e %+1.16e\n',xcoord(j),ycoord(j),Level);
    end
    i = i+N+1;
end
fclose(fid);

%% Done
filename = strcat('pseudospectrum_',strtrim(method),'.eps');
print(1,'-depsc2',filename);
