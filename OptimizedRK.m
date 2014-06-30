clear all;
clf;

fprintf('Make sure you have run cvx_setup!\n');
path(path,'/home/ghosh/Work/Codes/RK-opt/polyopt');
path(path,'/home/ghosh/Work/Codes/RK-opt/RKtools');
path(path,'/home/ghosh/Work/Codes/RK-opt/RK-coeff-opt');

filename=input('Enter eigenvalues filename: ','s');
x = load(filename);
z = x(:,1) + 1i*x(:,2);

p = input('Enter order of method : ');
s = input('Enter number of stages: ');

crit = input('Optimize for accuracy (acc) or SSP coeff (ssp): ','s');

[h,coeff]=opt_poly_bisect(z,s,p,'chebyshev','do_plot',true);
rk = rk_opt(s,p,'erk',crit,'poly_coeff_ind',((p+1):s), ...
    'poly_coeff_val',coeff((p+2):s+1),'display','iter');

plot(real(h*z),imag(h*z),'bo');
hold on;
plotstabreg_func(coeff,1);
hold off;

fprintf('Maximum CFL: %1.16E\n',h);
fprintf('A:\n');
fprintf([repmat('%1.16E\t', 1, size(rk.A, 2)) '\n'], rk.A');
fprintf('b:\n');
fprintf([repmat('%1.16E\t', 1, size(rk.b, 2)) '\n'], rk.b');
fprintf('c:\n');
fprintf([repmat('%1.16E\t', 1, size(rk.c, 2)) '\n'], rk.c');