clear all;
close all;

fprintf('Make sure you have run cvx_setup and set the path correctly!\n');

filename=input('Enter eigenvalues filename: ','s');
x = load(filename);
z = x(:,1) + 1i*x(:,2);

p    = input('Enter order of method : ');
smin = input('Enter minimum number of stages: ');
smax = input('Enter maximum number of stages: ');
crit = input('Optimize for accuracy (acc) or SSP coeff (ssp): ','s');

style1='-r';
style2='b.';
ns1=3;
ns2=ns1-1;

max_name_size=16;
legend_str=char(zeros(2,max_name_size));
count=0;

scrsz = get(0,'ScreenSize');
figID = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
for s=smin:smax
    [h,coeff]=opt_poly_bisect(z,s,p,'chebyshev','do_plot',false);
    rk = rk_opt(s,p,'erk',crit,'poly_coeff_ind',((p+1):s), ...
        'poly_coeff_val',coeff((p+2):s+1),'display','iter');
    if (~isstruct(rk))
        fprintf('Unable to find RK method of order %d with %d stages ', ...
            p,s);
        fprintf(' (%s optimized).\n',crit);
    else
        count = count+1;
        figure(figID);
        plot(real(h*z),imag(h*z),style2,'MarkerSize',8);
        hold on;
        [xa,ya,R] = StabilityRegion(coeff,1);
        contour(xa,ya,R,[1 1],style1);
        string1=['hl   (',sprintf('%2d',s),' stages)'];
        string2=['R(z) (',sprintf('%2d',s),' stages)'];
        legend_str(1,:)=string1;
        legend_str(2,:)=string2;
        xlabel('Real','FontName','Times','FontSize',14);
        ylabel('Imaginary','FontName','Times','FontSize',14);
        set(gca,'FontName','Times','FontSize',10);
        legend(legend_str(1:2*count,:),'Location','BestOutside');
        grid on;
        hold off;
        % save this figure
        figname = strcat('rk_opt_',sprintf('%1d',p),'_',sprintf('%02d',s));
        saveas(figID,strcat(figname,'.fig'),'fig');
        print(figID,'-depsc2',strcat(figname,'.eps'));
        
        % create a name for this method
        name = strcat('rk_opt_',sprintf('%1d',p),'_',sprintf('%02d',s));
        fprintf('Writing RK method to file %s.\n',strcat(name,'.inp'));
        fileID = fopen(strcat(name,'.inp'),'w');
        fprintf(fileID,'1\n');
        fprintf(fileID,'begin\n');
        fprintf(fileID,'\tname         %s\n',name);
        fprintf(fileID,'\tclass        rk\n');
        fprintf(fileID,'\tnstages      %d\n',s);
        fprintf(fileID,'\torder        %d\n',p);
        fprintf(fileID,'\thmax         %1.16e\n',h);
        fprintf(fileID,'\tA\n');
        for ii = 1:size(rk.A,1)
            fprintf(fileID,'\t\t%1.16E\t',rk.A(ii,:));
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'\tb\n');
        fprintf(fileID,'\t\t%1.16E\t',rk.b(:));
        fprintf(fileID,'\n');
        fprintf(fileID,'end\n');
        fclose(fileID);
    end
end

