%Creo matrici per script fortran
close all; clear; clc;

cr =1;
ct=1;
b=10;
m_piu1=31;
n_piu1=51;


Wing = build_wing(cr,ct,b,m_piu1,n_piu1, 2412, 0, 0);

surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
xlabel('chord')
ylabel('span')
%axis equal

x=Wing(:,:,1);
y=Wing(:,:,2);
z=Wing(:,:,3);

save X x
save Y y
save Z z

