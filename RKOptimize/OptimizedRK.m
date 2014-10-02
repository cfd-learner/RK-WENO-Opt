clear all;
close all;

fprintf('Make sure you have run cvx_setup and set the path correctly!\n');

filename=input('Enter eigenvalues filename: ','s');
x = load(filename);
z = x(:,1) + 1i*x(:,2);

p = input('Enter order of method : ');
smin = input('Enter minimum number of stages: ');
smax = input('Enter maximum number of stages: ');
crit = input('Optimize for accuracy (acc) or SSP coeff (ssp): ','s');

style1=['k: ','b: ','r: ','k-.','b-.','r-.'];
style2=['ko ','bo ','ro ','ks ','bs ','rs '];
ns1=3;
ns2=ns1-1;

max_name_size=15;
legend_str=char(zeros(2*(smax-smin+1),max_name_size));
count=0;
figure('units','normalized','position',[.1 .1 .6 .6])
for s=smin:smax
    [h,coeff]=opt_poly_bisect(z,s,p,'chebyshev','do_plot',false);
    rk = rk_opt(s,p,'erk',crit,'poly_coeff_ind',((p+1):s), ...
        'poly_coeff_val',coeff((p+2):s+1),'display','iter');
    
    ss=s-smin+1;
    plot(real(h*z),imag(h*z),strtrim(style2(ns1*ss-ns2:ns1*ss)), ...
        'MarkerSize',12);
    hold on;
    [xa,ya,R] = StabilityRegion(coeff,1);
    contour(xa,ya,R,[1 1],strtrim(style1(ns1*ss-ns2:ns1*ss)));
    hold on;
    string1=['hl   (',num2str(s,'%1d'),' stages)'];
    string2=['R(z) (',num2str(s,'%1d'),' stages)'];
    legend_str(2*ss-1,:)=string1;
    legend_str(2*ss,:)  =string2;
    
    % create a name for this method
    name = strcat('rk_opt_',num2str(p,'%1d'),'_',num2str(s,'%2d'));
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

xlabel('Real','FontName','Times','FontSize',15);
ylabel('Imaginary','FontName','Times','FontSize',15);
set(gca,'FontName','Times','FontSize',15);
legend(legend_str(1:2*(smax-smin+1),:),'Location','BestOutside');
grid on;
hold off;
