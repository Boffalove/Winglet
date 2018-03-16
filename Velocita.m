function [V, vv, Vh, Vv]=Velocita(Xa,X,Gamma)
%
% Questa funzione calcola la velocità indotta dalla griglia vorticosa Xa nel punto
% X. Gamma è la matrice i cui elementi sono le circolazioni degli anelli vorticosi
% della griglia.
% V contiene le velocità indotte da ciascun anello della griglia nel punto X
% vv è il vettore velocità indotta da tutta la griglia nel punto X
% Vh contiene le velocità indotte da ciascun segmento orizzontale della griglia
% (supposto a circolazione unitaria) nel punto X
% Vv è l'analogo di Vh ma per i segmenti verticali.

% [m,n,3]=size(Xa)
%
% V : (m-1, n-1)
% v : (3, 1)
% Vh: (m, n-1)
% Vv: (m-1,n)
% X : (3,1,3)
% Gamma: (m-1, n-1)
%
%===============================================================================

ma = size(Xa,1);
na = size(Xa,2);

vv = zeros(3,1);
Vv = zeros(ma-1, na,3);
Vh = zeros(ma, na-1,3);

for k = 1 : na-1
  
  Vv(:,k,:) = induced_line_speed([Xa(1:end-1,k,:), Xa(2:end,k,:)], ones(ma-1,1), X);
  Vh(:,k,:) = induced_line_speed([Xa(:,k+1,:), Xa(:,k,:)],ones(ma,1),X);
  
end

Vv(:,na,:) = induced_line_speed([Xa(1:end-1,na,:), Xa(2:end,na,:)], ones(ma-1,1), X);


[V,v]= Assemblaggio(Vv,Vh,Gamma);
vv(:)=v(1,1,:);


end





%%%%%Funzione Assemblaggio%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [V, vp] = Assemblaggio(Vv, Vh, Gamma)


m = size(Vv,1);
n = size(Vh,2);

V = zeros(m,n,3);

for k = 1 : n
  
   V(:,k,:) = Vv(:,k,:)-Vv(:,k+1,:) + Vh(1:end-1,k,:)-Vh(2:end,k,:);
   
end

V = V.*Gamma;

vp = sum(sum(V,1),2);

end

