clear all;
curr_folder = pwd;

fprintf('Initial script to set up CVX and add RK-opt paths\n');
cvx_path=input('Enter path to cvx: ', 's');
rk_opt_path=input('Enter path to RK Optimizer: ', 's');

cd(cvx_path);
cvx_setup;
 
cd(curr_folder);
str=strcat(rk_opt_path,'polyopt');
path(path,str);
str=strcat(rk_opt_path,'RKtools');
path(path,str);
str=strcat(rk_opt_path,'RK-coeff-opt');
path(path,str);
 