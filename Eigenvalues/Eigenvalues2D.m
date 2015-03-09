% Script to run HyPar to simulate the following:
% Case: Advection of a Sine Wave
% Model: 2D Linear Avection

clear all;
close all;

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');

% Compile the code to generate the initial and exact solutions
cmd = ['g++ ',hypar_path, ...
    '/Examples/2D/LinearAdvection/SineWave/aux/exact.c ', ...
    '-o EXACT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

% add the Matlab scripts directory in HyPar to path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Get the default
[~,~,~,~,~,~,~,~,~,~,~,~,par_scheme,~, ...
 ~,~,~,~, ~,input_mode, ...
 output_mode,n_io,~,~,~,~,~,~,~, ...
 ~,~,~,~,~,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
 verbose] = SetDefaults();
par_type = 'nonconservative-2stage';

% set problem specific input parameters
ndims = 2;
nvars = 1;
iproc = [1 1];
ghost = 3;

% set grid size;
N = [20 20];

% specify spatial discretization scheme
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
hyp_flux_split  = 'no';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 0;
nl      = 1;
eps     = 1e-6;

% time integration
dt      = 0.04;
t_final = 1.0;
niter   = int32(t_final/dt);

% set physical model and related parameters
model     = 'linear-advection-diffusion-reaction';
N_angles  = 100;
wave_angles = 0.0:2*pi/N_angles:2*pi;

% other options
cons_check      = 'yes';
screen_op_iter  = 1;
op_format       = 'none';
op_overwrite    = 'yes';
file_op_iter    = int32((t_final/1)/dt);
ip_type         = 'ascii';

% set time-integration scheme
ts      = 'rk';
tstype  = 'ssprk3';
use_petsc = 0;

% set boundaries
nb = 4;
bctype = ['periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'];
bcdim  = [0; 0; 1; 1;];
face   = [1; -1; 1; -1];
limits = [0 0 0.0 1.0; 0 0 0.0 1.0; 0.0 1.0 0 0; 0.0 1.0 0 0];

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
exact_exec = './EXACT > exact.log 2>&1';
hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ',hypar, ...
              ' > run.log 2>&1 '];

FFunction_eval_fname  = strcat('2D_',hyp_scheme,'.dat');
if (exist(FFunction_eval_fname,'file'))
    system(['rm ',FFunction_eval_fname]);
    fprintf('Removed existing file %s.\n',FFunction_eval_fname);
end
eval_fid = fopen(FFunction_eval_fname,'w');

for theta = wave_angles
    advection = [cos(theta) sin(theta)];
    fprintf('Advection wave angle: %4.2f\n',theta);
    % write the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts,tstype, ...
        hyp_scheme,hyp_flux_split,hyp_int_type,par_type,par_scheme, ...
        dt,cons_check,screen_op_iter,file_op_iter,op_format,ip_type, ...
        input_mode,output_mode,n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,bcdim,face,limits);
    WritePhysicsInp_LinearADR(advection);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    % Generate initial and exact solution
    system(exact_exec);
    % Run the simulation
    system(hypar_exec);

    % clean up
    system('rm -rf errors.dat conservation.dat function_counts.dat *.inp *.log *.bin *.eps');

    % Read operator matrix and compute it's eigenvalues
    fname_FFunction_root   = 'Mat_FFunction_';
    fname_extn = '.dat';

    index = sprintf('%05d',0);
    filename_FFunction = strcat(fname_FFunction_root,index,fname_extn);
    if (exist(filename_FFunction,'file'))
        fid = fopen(filename_FFunction,'r');
        ndof = fscanf(fid,'%d',1);
        fprintf('ndof = %d, ',ndof);
        A = zeros(ndof,ndof);
        nnz = 0;
        while (~feof(fid))
            coord = fscanf(fid,'%d',2);
            if (max(size(coord)) > 0)
                nnz = nnz + 1;
                A(coord(1),coord(2)) = fscanf(fid,'%f',1);
            end
        end
        fprintf('nnz = %d.\n',nnz);
        fclose(fid);
        FFunction_Mat = sparse(A);

        fprintf('  Computing %5d eigenvalues of FFunction matrix... ',ndof);
        tic;
        lambdaFFunction = eig(full(FFunction_Mat));
        waqt = toc;
        fprintf('%f seconds.\n',waqt);
        fprintf('  Saving eigenvalues to file %s.\n',FFunction_eval_fname);
        for n = 1:ndof
            fprintf(eval_fid,'%+1.16e %+1.16e\n', ...
                    real(lambdaFFunction(n:n)), ...
                    imag(lambdaFFunction(n:n)));
        end

        figure(1);
        plot(real(lambdaFFunction),imag(lambdaFFunction),'b.');
        hold on;
    else
        fprintf('File %s does not exist.\n',filename_FFunction);
    end
    system('rm -rf Mat*.dat');
end

fclose(eval_fid);
system('rm -rf EXACT');
