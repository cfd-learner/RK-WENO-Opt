clear all;
close all;

SystemName = '1D Euler';
CaseName = 'Sod Shock Tube';

% Ask for path to HyPar source directory
hypar_path = input('Enter path to HyPar source: ','s');
ti_path = input('Enter path to TI Methods: ','s');

% Ask for initial value of dt and tolerance
dt_init = input('Enter initial dt: ');
tolerance = input('Enter step size tolerance: ');

% Add to MATLAB path
path(path,strcat(hypar_path,'/Examples/Matlab/'));

% Compile the code to generate the initial solution
cmd = ['g++ ',hypar_path, ...
       '/Examples/1D/Euler1D/SodShockTube/aux/init.C ', ...
       '-o INIT'];
system(cmd);
% find the HyPar binary
hypar = [hypar_path,'/bin/HyPar'];

file_extn = '.inp';
orders = [2,3,4];

% Get the default
[~,~,~,~,~,~,~,~,~, ...
    hyp_flux_split,~,par_type,par_scheme,~,cons_check, ...
    screen_op_iter,file_op_iter,~,~,input_mode, ...
    output_mode,n_io,op_overwrite,~,~,~,~,~,~, ...
    mapped,borges,yc,nl,eps,p,rc,xi,wtol,lutype,norm,maxiter,atol,rtol, ...
    verbose] = SetDefaults();

% Set initial solution file type to ASCII
ip_type = 'ascii';

% set problem specific input parameters
ndims = 1;
nvars = 3;
iproc = 1;
ghost = 3;

% use TI methods optimized for SSP or accuracy?
opt_type = 'SSP';

% specify a nice, high-order spatial discretization scheme
hyp_scheme = 'weno5';
hyp_int_type = 'characteristic';

% set final time
t_final = 0.2;

% set physical model and related parameters
model = 'euler1d';
gamma = 1.4;
grav  = 0.0;
upw   = 'roe';

% set boundaries
nb = 2;
bctype = ['extrapolate'; 'extrapolate'];
bcdim = [0; 0];
face = [1; -1];
limits = [0 0; 0 0];

% plotting styles
colors = ['k','b','r'];
pointstyle = ['o','s','d','^','v','<','>','p','h','x','+','*'];
linestyle = '--';

% turn off solution output to file
op_format = 'text';

% set grid size;
N = 201;

% set the commands to run the executables
init_exec = './INIT > init.log 2>&1';
clean_exec = 'rm -rf *.inp *.dat *.log';

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
ti_path = [ti_path,sprintf('%1d',ndims),'D'];
if (strcmp(hyp_scheme,'weno5'))
    ti_path = [ti_path,'_','NonCompact5'];
elseif (strcmp(hyp_scheme,'crweno5'))
    ti_path = [ti_path,'_','Compact5'];
else
    fprintf('Incorrect hyp_scheme specified.\n');
    return;
end
ti_path = [ti_path,'_NonLinear_2S'];
count = 1;
n_o = 1;
for order = orders
    nstages = order:order+11;
    n_s = 0;
    TSStages = zeros(size(nstages,2),1);
    TSMaxDt  = zeros(size(nstages,2),1);
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
        fprintf('\tdt=%1.6e, dt/s=%1.6e, factor=%8.6f: ',dt,dt/stages,dt_factor);
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
        WritePhysicsInp_Euler1D(gamma,grav,upw);
        WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
        WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
        if (flag_opt)
            % create a sym link to the file containing the TI method
            system(['ln -sf ',filename,' time_method.inp']);
        end
        system(init_exec);
        system('ln -sf initial.inp exact.inp');
        hypar_exec = [hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                      ' > run.log 2>&1'];
        system(hypar_exec);
        [err,wt] = ReadErrorDat(ndims);
        [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
        % save the output file
        op_fname = strcat('op_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                          sprintf('%1d',order),'_', ...
                          sprintf('%02d',stages),'_', ...
                          'dt',sprintf('%10.8f',dt),'.txt');
        system(['mv op.dat ',op_fname]);
        log_fname = strcat('log_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                           sprintf('%1d',order),'_', ...
                           sprintf('%02d',stages),'_', ...
                           'dt',sprintf('%10.8f',dt),'.txt');
        system(['mv run.log ',log_fname]);
        system(clean_exec);
        if (~min(isfinite(err)))
            fprintf('failed %1.16e.\n',err(2));
            continue;
        end
        fprintf('passed %1.16e.\n',err(2));
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

            fprintf('\tdt=%1.6e, dt/s=%1.6e, factor=%8.6f: ',dt,dt_new/stages,dt_factor);
            niter = floor(t_final/dt_new);
            petscdt = [' -ts_dt ',num2str(dt_new,'%1.16e'),' '];
            petscft = [' -ts_final_time ',num2str((niter*dt_new),'%f'),' '];
            petscms = [' -ts_max_steps ',num2str(100*niter,'%d'),' '];
            WriteSolverInp(ndims,nvars,N,iproc,ghost,niter,strtrim(ts), ...
                strtrim(tstype),hyp_scheme,hyp_flux_split,hyp_int_type, ...
                par_type,par_scheme, dt_new,cons_check,screen_op_iter, ...
                file_op_iter,op_format,ip_type,input_mode,output_mode, ...
                n_io,op_overwrite,model);
            WriteBoundaryInp(nb,bctype,bcdim,face,limits);
            WritePhysicsInp_Euler1D(gamma,grav,upw);
            WriteWenoInp(mapped,borges,yc,nl,eps,p,rc,xi,wtol);
            WriteLusolverInp(lutype,norm,maxiter,atol,rtol,verbose);
            if (flag_opt)
                % create a sym link to the file containing the TI method
                system(['ln -sf ',filename,' time_method.inp']);
            end
            system(init_exec);
            system('ln -sf initial.inp exact.inp');
            hypar_exec = [hypar,' ',petsc_flags,petscdt,petscft,petscms, ...
                          ' > run.log 2>&1'];
            system(hypar_exec);
            [err,wt] = ReadErrorDat(ndims);
            [~,fcounts,~,~,~,~,~,~] = ReadFunctionCounts();
            % save the output file and log
            op_fname = strcat('op_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                              sprintf('%1d',order),'_', ...
                              sprintf('%02d',stages),'_', ...
                              'dt',sprintf('%10.8f',dt_new),'.txt');
            system(['mv op.dat ',op_fname]);
            log_fname = strcat('log_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                              sprintf('%1d',order),'_', ...
                              sprintf('%02d',stages),'_', ...
                              'dt',sprintf('%10.8f',dt_new),'.txt');
            system(['mv run.log ',log_fname]);
            system(clean_exec);
            if (~min(isfinite(err)))
                fprintf('failed. %1.16e.\n',err(2));
                dt_max = dt_new;
                dt_factor = 0.5 * dt_factor;
            else
                fprintf('passed %1.16e.\n',err(2));
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
    
        % set legend string
        name_str = [sprintf('%1d',order),'(',sprintf('%2d',stages),')'];
        legend_str(count,:) = name_str;

        count = count+1;
        n_s = n_s + 1;
    end
    
    TSMaxDt = TSMaxDt(1:n_s);
    TSStages = TSStages(1:n_s);
    % write to file
    data_fname = ['dtmax_',strtrim(hyp_scheme),'_',opt_type,'_', ...
                   sprintf('%1d',order),'.txt'];
    fprintf('Saving data to %s.\n',data_fname);
    data_fid = fopen(data_fname,'w');
    for ii = 1:n_s
        fprintf(data_fid,'%1.16e %1.16e\n', ...
                TSStages(ii),TSMaxDt(ii));
    end
    fclose(data_fid);
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

figname = [opt_type,'_',hyp_scheme,'_',sprintf('%1d',orders)];
if (count > 1)

    title_str = [SystemName,'(',CaseName,'): N=',sprintf('%d',N),', ', ...
                 'Spatial scheme ',hyp_scheme,', ', ...
                 'Opt type ',opt_type];

    figure(figErrDt);
    xlabel('dt','FontName','Times','FontSize',14,'FontWeight','normal');
    ylabel('Error','FontName','Times','FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',10,'FontName','Times');
    legend(legend_str,'Location','eastoutside');
    axis auto;
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
    axis auto;
    title(title_str);
    grid on;
    hold off;
    
    figure(figMaxDtStages);
    xlabel('Number of Stages','FontName','Times','FontSize',14,'FontWeight','normal');
    ylabel('Maximum stage time step','FontName','Times','FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',10,'FontName','Times');
    legend(legend_ptr,'Location','best');
    axis auto;
    title(title_str);
    grid on;
    hold off;

    % print figures to file
    print(figErrDt,'-depsc2',['fig_',figname,'_ErrDt.eps']);
    saveas(figErrDt,['fig_',figname,'_ErrDt.fig']);
    print(figErrCost,'-depsc2',['fig_',figname,'_ErrCost.eps']);
    saveas(figErrCost,['fig_',figname,'_ErrCost.fig']);
    print(figMaxDtStages,'-depsc2',['fig_',figname,'_MaxDtStages.eps']);
    saveas(figMaxDtStages,['fig_',figname,'_MaxDtStages.fig']);
end

% clean up
system('rm -rf INIT');
