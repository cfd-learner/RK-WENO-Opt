clear all;
close all;

SystemName = '2D Navier-Stokes';
CaseName = 'Isentropic Vortex Convection';

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');
ti_path = input('Enter path to TI Methods: ','s');

% Ask for initial value of dt and tolerance
dt_init = input('Enter initial dt: ');
tolerance = input('Enter step size tolerance: ');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['gcc ',hypar_path, ...
       '/Examples/2D/NavierStokes2D/InviscidVortexConvection/aux/init.c ', ...
       '-lm -o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

file_extn = '.inp';
orders = [2,3,4];
% use TI methods optimized for SSP or accuracy?
opt_type = 'SSP';

% Get the default
[~,~,~,~,~,~,~,~,~, ...
    ~,~,~,par_scheme,~,~, ...
    ~,~,~, ~,input_mode, ...
    output_mode,n_io,~,~,~,~,~,~,~, ...
    ~,~,~,~,~,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% set problem specific input parameters
ndims = 2;                 % number of space dimensions
nvars = 4;                 % number of variables in state vector
iproc = [2 2];             % number of processors in each dimension
ghost = 3;                 % number of ghost points
N = [32 32];               % grid dimensions

% final time
t_final = 20.0;

% specify spatial discretization scheme details
hyp_scheme      = 'weno5';
hyp_int_type    = 'components';
hyp_flux_split  = 'no';
par_type        = 'nonconservative-2stage';
% parameters controlling the WENO-type schemes
mapped  = 0;
borges  = 0;
yc      = 0;
nl      = 1;
eps     = 1e-6;

% set physical model and related parameters
model     = 'navierstokes2d';
gamma     = 1.4;                  % specific heat ratio
upw       = 'roe';                % choice of upwinding scheme
Prandtl   = 0.72;                 % Prandtl number
Reynolds  = -1.0;                 % Inviscid flow
Minf      = 1.0;                  % reference Mach number
grav      = [0.0 0.0];            % gravitational force vector
rho_ref   = 1.0;                  % reference altitude density
p_ref     = 1.0;                  % reference altitude pressure
HB        = 0;                    % type of hydrostatic balance
BV        = 0.0;                  % Brunt-Vaisala frequency
GasConst  = 287.058;              % Universal gas constant
% other options
cons_check      = 'yes';
screen_op_iter  = 10;
op_overwrite    = 'yes';
file_op_iter    = 999999;
ip_type         = 'binary';
op_format       = 'none';

% set boundaries
nb = 4;
bctype = ['periodic'; ...
          'periodic'; ...
          'periodic'; ...
          'periodic'];
bcdim     = [0; 0; 1; 1;];
face      = [1; -1; 1; -1];
limits    = [0 0 0 10.0; 0 0 0 10.0; 0 10.0 0 0; 0 10.0 0 0];

maxerr = 1.0;

% plotting styles
colors = ['k','b','r'];
pointstyle = ['o','s','d','^','v','<','>','p','h','x','+','*'];
linestyle = '--';

% set the commands to run the executables
nproc = 1;
for i = 1:max(size(iproc))
    nproc = nproc * iproc(i);
end
init_exec  = './INIT  > init.log  2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log op.*';

% Generate or find the reference solution
RefFlag = input('Generate reference solution? ','s');
if (strcmp(RefFlag,'yes'))
    fprintf('Generating reference solution...\n');
    % small time step for reference solution
    dt_ref = 0.05 * dt_init;    
    % use explicit RK 4-stage, 4th order
    ts_ref = 'rk';              
    tstype_ref = '44';          
    niter = int32(t_final/dt_ref);
    % set reference solution output type as binary
    op_format = 'binary';
    % Write out the input files for HyPar
    WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,ts_ref, ...
        tstype_ref,hyp_scheme,hyp_flux_split,hyp_int_type, ...
        par_type,par_scheme, dt_ref,cons_check,screen_op_iter, ...
        file_op_iter,op_format,ip_type,input_mode,output_mode, ...
        n_io,op_overwrite,model);
    WriteBoundaryInp(nb,bctype,bcdim,face,limits);
    WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                   rho_ref,p_ref,HB,BV,GasConst);
    WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
    WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
    % Generate the initial solution
    system(init_exec);
    % Run HyPar
    hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                  hypar,' > run.log 2>&1'];
    system(hypar_exec);
    % convert the output to the input format
    system(['gcc ',hypar_path,'/Extras/BinaryOPToInitialSolution.c ', ...
            '-o BINOP2INP']);
    fid = fopen('bin.inp','w');
    fprintf(fid,'op.bin');
    fclose(fid);
    system('./BINOP2INP < bin.inp 2>&1 > conv.log && rm bin.inp');
    system('mv solution.inp reference.bin');
    system('mv run.log reference.log');
    % save the reference solution and log in a separate directory
    dir_name = strcat('refsoln_',sprintf('%03d_',N),hyp_scheme, ...
                      '_',sprintf('%04.1f',t_final));
    if (exist(dir_name,'file'))
        fprintf('Removing existing directory %s.\n',dir_name);
        system(['rm -rf ',dir_name]);
    end
    mkdir(dir_name);
    system(['mv op.bin *.inp reference.bin reference.log ',dir_name]);
    system(['rm -rf ',dir_name,'/initial.inp']);
    % clean up
    system('rm -rf *.inp *.log *.dat *.bin BINOP2INP');
    % create a shortcut in the current folder for the reference solution
    system(['ln -sf ',dir_name,'/reference.bin reference.bin']);
else
    refpath = input('Enter path to reference solution: ','s');
    if (~strcmp(refpath,'./')) && (~strcmp(refpath,'.'))
        if (~exist([refpath,'/reference.bin'],'file'))
            fprintf('Error: reference solution not found in %s.\n', ...
                    refpath);
            return;
        end
        system(['ln -sf ',refpath,'/reference.bin reference.bin']);
    end
end

% open figure window
scrsz = get(0,'ScreenSize');
figErrDt   = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
figErrCost = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
figMaxDtStages = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% initialize legend string
legend_str = char(zeros(size(orders,2)*20,5));
legend_ptr = char(zeros(size(orders,2),9));

% set maximum number of data points
ref_levels = 1000;

ti_path = [ti_path,'/',opt_type,'-optimized/'];
ti_path = [ti_path,sprintf('%1d',1),'D'];
if (strcmp(hyp_scheme,'weno5'))
    ti_path = [ti_path,'_','NonCompact5'];
elseif (strcmp(hyp_scheme,'crweno5'))
    ti_path = [ti_path,'_','Compact5'];
else
    fprintf('Incorrect hyp_scheme specified.\n');
    return;
end
count = 1;
MinDt   = zeros(size(orders,2)*20,1);
MaxDt   = zeros(size(orders,2)*20,1);
MinErr  = zeros(size(orders,2)*20,1);
MaxErr  = zeros(size(orders,2)*20,1);
MinCost = zeros(size(orders,2)*20,1);
MaxCost = zeros(size(orders,2)*20,1);
n_o = 1;
for order = orders
    nstages = order:order+11;
    %if (strcmp(hyp_scheme,'crweno5'))
    %    if (order == 2)
    %        nstages = [2,4];
    %    elseif (order == 3)
    %        nstages = [3,6];
    %    elseif (order == 4)
    %        nstages = [4,9];
    %    end
    %elseif (strcmp(hyp_scheme,'weno5'))
    %    if (order == 2)
    %        nstages = [2,5];
    %    elseif (order == 3)
    %        nstages = [3,6];
    %    elseif (order == 4)
    %        nstages = [4,9];
    %    end
    %end
    n_s = 0;
    TSStages = zeros(size(nstages,2),1);
    TSMaxDt = zeros(size(nstages,2),1);
    for stages = nstages
        ts = 'rk';
        if (stages == order)
            if (order == 2)
                tstype = '2a';
            elseif (order == 3)
                tstype = '3';
            elseif (order == 4)
                tstype = '4';
            end
            flag_opt = 0;
        else
            tstype = ['rk_opt_',sprintf('%1d',order),'_',sprintf('%02d',stages)];
            flag_opt = 1;
        end
        fprintf('Order %2d, Stages %2d, TS %s, TSType %s:\n', ...
                 order,stages,ts,tstype);
        % set PETSc time-integration flags (comment to turn off)
        petsc_flags = sprintf('%s', ...
            '-use-petscts ', ...
            '-ts_type ',strtrim(ts),' ', ...
            '-ts_',strtrim(ts),'_type ',strtrim(tstype),' ', ...
            '-ts_adapt_type none ', ...
            ' ');
        
        % check if this method exists
        if (flag_opt)
            filename = [ti_path,'/',tstype,file_extn];
            if (~exist(filename,'file'))
                fprintf('not found %s.\n',filename);
                continue;
            end
        end
    
        % set dt to its initial value
        dt = dt_init;
        dt_factor = 1.0;
        r = 1;

        % set max time step size to final time
        dt_max = t_final;
    
        % preallocate arrays for dt, error, wall times and function counts
        TimeStep  = zeros(ref_levels,1);
        Errors    = zeros(ref_levels,3);
        Walltimes = zeros(ref_levels,2);
        FCounts   = zeros(ref_levels,1);

        % run simulation with initial dt
        fprintf('\tdt=%1.6e, factor=%8.6f: ',dt,dt_factor);
        niter = floor(t_final/dt);
                petscdt = [' -ts_dt ',num2str(dt,'%1.16e'),' '];
        petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
        petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
        WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts), ...
                strtrim(tstype),hyp_scheme,hyp_flux_split,hyp_int_type, ...
                par_type,par_scheme, dt,cons_check,screen_op_iter, ...
                file_op_iter,op_format,ip_type,input_mode,output_mode, ...
                n_io,op_overwrite,model);
        WriteBoundaryInp(nb,bctype,bcdim,face,limits);
        WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                       rho_ref,p_ref,HB,BV,GasConst);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        if (flag_opt)
            % create a sym link to the file containing the TI method
            system(['ln -sf ',filename,' time_method.inp']);
        end
        system(init_exec);
        system('ln -sf reference.bin exact.inp');
        hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                      hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1' ...
                     ];
        system(hypar_exec);
        [err,wt] = ReadErrorDat(ndims);
        [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
        system(clean_exec);
        if (~min(isfinite(err)))
            fprintf('failed.\n');
            continue;
        end
        fprintf('passed. err=%1.4e\n',err(2));
        TimeStep(r)     = dt;
        Errors(r,:)     = err;
        Walltimes(r,:)  = wt;
        FCounts(r)      = fcounts;
        r = r+1;

        while ((dt_factor > tolerance) && (r < ref_levels))
            dt_new = dt * (1.0+dt_factor);
            while (dt_new > dt_max)
                dt_factor = 0.5 * dt_factor;
                dt_new = dt * (1.0+dt_factor);
            end

            % estimate error from previous error based on theoretical order
            err_theoretical = Errors(r-1,2) * (dt_new/dt)^order;

            fprintf('\tdt=%1.6e, factor=%8.6f: ',dt_new,dt_factor);
            niter = floor(t_final/dt_new);
            petscdt = [' -ts_dt ',num2str(dt_new,'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str(t_final,'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
            WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts), ...
                strtrim(tstype),hyp_scheme,hyp_flux_split,hyp_int_type, ...
                par_type,par_scheme, dt_new,cons_check,screen_op_iter, ...
                file_op_iter,op_format,ip_type,input_mode,output_mode, ...
                n_io,op_overwrite,model);
            WriteBoundaryInp(nb,bctype,bcdim,face,limits);
            WritePhysicsInp_NavierStokes2D(gamma,upw,Prandtl,Reynolds,Minf,grav, ...
                                           rho_ref,p_ref,HB,BV,GasConst);
            WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
            WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
            if (flag_opt)
                % create a sym link to the file containing the TI method
                system(['ln -sf ',filename,' time_method.inp']);
            end
            system(init_exec);
            system('ln -sf reference.bin exact.inp');
            hypar_exec = ['$MPI_DIR/bin/mpiexec -n ',num2str(nproc),' ', ...
                          hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1' ...
                          ];
            system(hypar_exec);
            [err,wt] = ReadErrorDat(ndims);
            [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
            system(clean_exec);
            if ((~min(isfinite(err))) || (err(2)/err_theoretical > 2))
                fprintf('failed.\n');
                dt_max = dt_new;
                dt_factor = 0.5 * dt_factor;
            else
                fprintf('passed. err=%1.4e\n',err(2));
                dt = dt_new;
                TimeStep(r) = dt;
                Errors(r,:) = err;
                Walltimes(r,:) = wt;
                FCounts(r) = fcounts;
                r = r+1;
            end
        end
    
        % Isolate the L2 Error
        L2Errors = Errors(:,2);
        
        TSStages(n_s+1) = stages;
        TSMaxDt(n_s+1) = max(TimeStep(1:r-1))/stages;
    
        % To be used in setting axis limits
        MinDt(count)   = min(TimeStep(1:r-1));
        MaxDt(count)   = max(TimeStep(1:r-1));
        MinErr(count)  = min(L2Errors(1:r-1));
        MaxErr(count)  = max(L2Errors(1:r-1));
        MinCost(count) = min(FCounts (1:r-1));
        MaxCost(count) = max(FCounts (1:r-1));

        if (flag_opt && n_s)
            style = [linestyle,colors(n_o),pointstyle(n_s)];
            linewidth = 1;
        else
            style = ['-',colors(n_o)];
            linewidth = 2;
        end
        % plot errors
        figure(figErrDt);
        loglog(TimeStep(1:r-1),L2Errors(1:r-1),style,'linewidth',linewidth, ...
            'MarkerSize',8);
        hold on;
        % plot cost
        figure(figErrCost);
        loglog(FCounts(1:r-1),L2Errors(1:r-1),style,'linewidth',linewidth, ...
               'MarkerSize',8);
        hold on;
    
        % write to file
        data_fname = ['data_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                       sprintf('%1d',order),'_', ...
                       sprintf('%02d',stages),'.txt'];
        fprintf('Saving data to %s.\n',data_fname);
        data_fid = fopen(data_fname,'w');
        for ii = 1:r-1
            fprintf(data_fid,'%1.16e %1.16e %1.16e\n', ...
                    TimeStep(ii), ...
                    L2Errors(ii), ...
                    FCounts(ii));
        end
        fclose(data_fid);
    
        % set legend string
        name_str = [sprintf('%1d',order),'(',sprintf('%2d',stages),')'];
        legend_str(count,:) = name_str;
    
        count = count+1;
        n_s = n_s + 1;
    end
    
    TSMaxDt = TSMaxDt(1:n_s);
    TSStages = TSStages(1:n_s);
    % plot max dt as a function of number of stages
    style = ['-',colors(n_o),'o'];
    figure(figMaxDtStages);
    plot(TSStages,TSMaxDt,style,'linewidth',1,'MarkerSize',8);
    if (order == 1)
        blah = 'st';
    elseif (order == 2)
        blah = 'nd';
    elseif (order == 3)
        blah = 'rd';
    else
        blah = 'th';
    end
    legend_ptr(n_o,:) = [sprintf('%1d',order),blah,' Order'];
    hold on;
    
    n_o = n_o + 1;
end

MinDt   = MinDt(1:count-1);
MaxDt   = MaxDt(1:count-1);
MinErr  = MinErr(1:count-1);
MaxErr  = MaxErr(1:count-1);
MinCost = MinCost(1:count-1);
MaxCost = MaxCost(1:count-1);
legend_str = legend_str(1:count-1,:);

figname = [opt_type,'_',hyp_scheme,'_',sprintf('%1d',orders)];
if (count > 1)

    title_str = [SystemName,'(',CaseName,'): N=',sprintf('%d^2',N(1)),', ', ...
                 'Spatial scheme ',hyp_scheme,', ', ...
                 'Opt type ',opt_type];

    figure(figErrDt);
    xlabel('dt','FontName','Times','FontSize',14,'FontWeight','normal');
    ylabel('Error','FontName','Times','FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',10,'FontName','Times');
    legend(legend_str,'Location','eastoutside');
    axis([min(MinDt)/2 max(MaxDt)*2 min(MinErr)/2 min(2*max(MaxErr),maxerr)]);
    title(title_str);
    grid on;
    hold off;

    figure(figErrCost);
    ylabel('Error (L_2)','FontName','Times','FontSize',14, ...
           'FontWeight','normal');
    xlabel('Number of RHS function calls','FontName','Times', ...
           'FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',10,'FontName','Times');
    legend(legend_str,'Location','eastoutside');
    axis([min(MinCost)/2 2*max(MaxCost) min(MinErr)/2 min(2*max(MaxErr),maxerr)]);
    title(title_str);
    grid on;
    hold off;
    
    figure(figMaxDtStages);
    xlabel('Number of Stages','FontName','Times','FontSize',14,'FontWeight','normal');
    ylabel('Maximum stage time step','FontName','Times','FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',10,'FontName','Times');
    legend(legend_ptr,'Location','best');
    title(title_str);
    grid on;
    hold off;

    % print figures to file
    print(figErrDt,'-depsc2',['fig_',figname,'_ErrDt.eps']);
    print(figErrCost,'-depsc2',['fig_',figname,'_ErrCost.eps']);
    print(figMaxDtStages,'-depsc2',['fig_',figname,'_MaxDtStages.eps']);
    
    % save figures
    saveas(figErrDt,['fig_',figname,'_ErrDt.fig'],'fig');
    saveas(figErrCost,['fig_',figname,'_ErrCost.fig'],'fig');
    saveas(figMaxDtStages,['fig_',figname,'_MaxDtStages.fig'],'fig');
end

% clean up
system('rm -rf INIT reference.bin');
