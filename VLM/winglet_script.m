%Creo matrici per script fortran
close all; clear; clc;

%dati
cr = 1;
ct = 1;
b = 10;

m_piu1 = 31;
n_piu1 = 11;

m = m_piu1-1;
n = n_piu1-1;

alpha = 5;
Re = 1e+6;

% costruzione ala
Wing = build_wing(cr,ct,b,m_piu1,n_piu1, 2412, 0, 0);

% surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
% xlabel('chord')
% ylabel('span')
%axis equal

x=Wing(:,:,1);
y=Wing(:,:,2);
z=Wing(:,:,3);

% file delle coordinate
save X.dat x -ascii
save Y.dat y -ascii
save Z.dat z -ascii

% file data.dat
id = fopen('Data.dat','W');
form = '%s=%d\n%s=%d\n%s=%f\n%s=%f\n%s=%d\n';
fprintf(id,form,'m',m,'n',n,'b',b,'alpha',alpha,'Re',Re);
fclose(id)

if ~exist('./subroutines.mod','file')
    mainid = fopen('Main.f90');
    text = textread('Main.f90','%s\t','whitespace','','delimiter','\n');
    fclose(mainid);
    
    mainid = fopen('Main.f90','W');
    
    idx = strcmp(text,'MODULE Subroutines');
    idx = find(idx==1);
    
    module = text(idx(1):end);
    program = text(1:idx(1)-1);
    
    text_new =[module; program];
    
    fprintf(mainid,'%s\n',text_new{:});
    fclose(mainid);
    
    system('gfortran -c Main.f90');
    mainid = fopen('Main.f90','W');
    
    fprintf(mainid,'%s\n',text{:});
    fclose(mainid);
end

system('gfortran -c Main.f90')
system('gfortran -o Main Main.f90 -L/usr/local/lib -llapack -lblas')
system('./Main')

