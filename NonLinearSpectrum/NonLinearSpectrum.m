clear all;
close all;

%% CRWENO5
data = load('wavenumber_CRWENO5JS.dat');
eigenvalue = zeros(9*size(data,1),1);
for i=1:size(data,1)
    ax = data(i,5);
    ay = data(i,4);
    eigenvalue(9*i+1) = data(i,3) + 1i*data(i,2);
    eigenvalue(9*i+2) = data(i,3) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+3) = data(i,3) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+4) = (data(i,3)+ax) + 1i*data(i,2);
    eigenvalue(9*i+5) = (data(i,3)+ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+6) = (data(i,3)+ax) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+7) = (data(i,3)-ax) + 1i*data(i,2);
    eigenvalue(9*i+8) = (data(i,3)-ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+9) = (data(i,3)-ax) + 1i*(data(i,2)-ay);
    
end
for i=1:size(eigenvalue,1)
    if (real(eigenvalue(i)) > 0)
        eigenvalue(i) = 0 + 1i*imag(eigenvalue(i));
    end
end
figure(1);
plot(eigenvalue,'ko');
grid on;
axis equal;
fid = fopen('eigenvalues_crweno5_1s.dat','w');
for i=1:size(eigenvalue,1)
    fprintf(fid,'%+1.16e %+1.16e\n',real(eigenvalue(i)),imag(eigenvalue(i)));
end
fclose(fid);

%% WENO5
data = load('wavenumber_WENO5JS.dat');
eigenvalue = zeros(9*size(data,1),1);
for i=1:size(data,1)
    ax = data(i,5);
    ay = data(i,4);
    eigenvalue(9*i+1) = data(i,3) + 1i*data(i,2);
    eigenvalue(9*i+2) = data(i,3) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+3) = data(i,3) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+4) = (data(i,3)+ax) + 1i*data(i,2);
    eigenvalue(9*i+5) = (data(i,3)+ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+6) = (data(i,3)+ax) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+7) = (data(i,3)-ax) + 1i*data(i,2);
    eigenvalue(9*i+8) = (data(i,3)-ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+9) = (data(i,3)-ax) + 1i*(data(i,2)-ay);
    
end
for i=1:size(eigenvalue,1)
    if (real(eigenvalue(i)) > 0)
        eigenvalue(i) = 0 + 1i*imag(eigenvalue(i));
    end
end
figure(2);
plot(eigenvalue,'ko');
grid on;
axis equal;
fid = fopen('eigenvalues_weno5_1s.dat','w');
for i=1:size(eigenvalue,1)
    fprintf(fid,'%+1.16e %+1.16e\n',real(eigenvalue(i)),imag(eigenvalue(i)));
end
fclose(fid);

%% CRWENO5
data = load('wavenumber_CRWENO5JS.dat');
eigenvalue = zeros(9*size(data,1),1);
for i=1:size(data,1)
    ax = 2*data(i,5);
    ay = 2*data(i,4);
    eigenvalue(9*i+1) = data(i,3) + 1i*data(i,2);
    eigenvalue(9*i+2) = data(i,3) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+3) = data(i,3) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+4) = (data(i,3)+ax) + 1i*data(i,2);
    eigenvalue(9*i+5) = (data(i,3)+ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+6) = (data(i,3)+ax) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+7) = (data(i,3)-ax) + 1i*data(i,2);
    eigenvalue(9*i+8) = (data(i,3)-ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+9) = (data(i,3)-ax) + 1i*(data(i,2)-ay);
    
end
for i=1:size(eigenvalue,1)
    if (real(eigenvalue(i)) > 0)
        eigenvalue(i) = 0 + 1i*imag(eigenvalue(i));
    end
end
figure(3);
plot(eigenvalue,'ko');
grid on;
axis equal;
fid = fopen('eigenvalues_crweno5_2s.dat','w');
for i=1:size(eigenvalue,1)
    fprintf(fid,'%+1.16e %+1.16e\n',real(eigenvalue(i)),imag(eigenvalue(i)));
end
fclose(fid);

%% WENO5
data = load('wavenumber_WENO5JS.dat');
eigenvalue = zeros(9*size(data,1),1);
for i=1:size(data,1)
    ax = 2*data(i,5);
    ay = 2*data(i,4);
    eigenvalue(9*i+1) = data(i,3) + 1i*data(i,2);
    eigenvalue(9*i+2) = data(i,3) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+3) = data(i,3) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+4) = (data(i,3)+ax) + 1i*data(i,2);
    eigenvalue(9*i+5) = (data(i,3)+ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+6) = (data(i,3)+ax) + 1i*(data(i,2)-ay);
    eigenvalue(9*i+7) = (data(i,3)-ax) + 1i*data(i,2);
    eigenvalue(9*i+8) = (data(i,3)-ax) + 1i*(data(i,2)+ay);
    eigenvalue(9*i+9) = (data(i,3)-ax) + 1i*(data(i,2)-ay);
    
end
for i=1:size(eigenvalue,1)
    if (real(eigenvalue(i)) > 0)
        eigenvalue(i) = 0 + 1i*imag(eigenvalue(i));
    end
end
figure(4);
plot(eigenvalue,'ko');
grid on;
axis equal;
fid = fopen('eigenvalues_weno5_2s.dat','w');
for i=1:size(eigenvalue,1)
    fprintf(fid,'%+1.16e %+1.16e\n',real(eigenvalue(i)),imag(eigenvalue(i)));
end
fclose(fid);