%Creo matrici per script fortran
close all; clear; clc;

%dati
cr = 1;
ct = 1;
b = 10;

m_piu1 = 11;
n_piu1 = 101;

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
save X.dat x -ascii;
save Y.dat y -ascii;
save Z.dat z -ascii;

% file data.dat
id = fopen('Data.dat','W');
form = '%s=%d\n%s=%d\n%s=%f\n%s=%f\n%s=%d\n';
fprintf(id,form,'m',m,'n',n,'b',b,'alpha',alpha,'Re',Re);
fclose(id);
clear id form
% ottengo la data di ultima modifica del file
main_f90_date = dir('Main.f90'); main_f90_date = main_f90_date.datenum;

if exist('./subroutines.mod','file')
    sub_date = dir('subroutines.mod'); sub_date = sub_date.datenum;
else
    sub_date = 0;
end

if exist('Main','file')
    main_date = dir('Main'); main_date = main_date.datenum;
else
    main_date = 0;
end
% se il file subroutines.mod non esiste oppure la data di ultima modifica è
% precedente a  quella di Main.f90 lo (ri)compilo

if ~exist('./subroutines.mod','file') || sub_date <= main_f90_date
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
    
    comp_sub = system('gfortran -c Main.f90');
    mainid = fopen('Main.f90','W');
    
    fprintf(mainid,'%s\n',text{:});
    fclose(mainid);
    clear mainid idx module program text text_new comp_sub
end

% se il file Main non esiste oppure la data di ultima modifica è
% precedente a  quella di Main.f90 lo (ri)compilo e lo linko

if (~exist('Main','file') || main_date <= main_f90_date) && exist('./subroutines.mod','file')
   comp_main = system('gfortran -c Main.f90');
   lnk_main = system('gfortran -o Main Main.f90 -L/usr/local/lib -llapack -lblas');
else
    comp_main = 0;
    lnk_main = 0;
end

clear main_date main_f90_date sub_date

% se Main compilato e linkato allora posso eseguire
if comp_main == 0 && lnk_main == 0
    % se si vuole l'output a schermo rimuovere 'eclispe >/dev/null
    execution = system('./Main eclipse >/dev/null');
else
    error('Main non compilato o non linkato')
end

clear lnk_main comp_main

if execution ~= 0
    error('********************************************* Failed execution *********************************************')
end

clear execution

% carico i file per il postprocess

if exist('Xc.dat','file');    load Xc.dat;    end
if exist('gamma.dat','file'); load gamma.dat; end
if exist('Vel.dat','file');   load Vel.dat;   end
if exist('Cl.dat','file');    load Cl.dat;    end
if exist('Cd.dat','file');    load Cd.dat;    end


