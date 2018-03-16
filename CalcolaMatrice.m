function A = CalcolaMatrice(Xa,Xc,N)

% Xa: griglia vorticosa.   Xa: (ma,na,3)
% Xc: Griglia punti di controllo.    Xc: (ma-1, na-1, 3)
% N:  Array contenente i versori normali ai pannelli.   N: (ma-1, na-1, 3)
% A:  Matrice dei coefficienti di influenza

% A(i,j) è la velocità (a circolazione unitaria) indotta dal vortice j nel punto
% i (*Ni)

m = size(Xc,1)
n = size(Xc,2)
Nc = m*n;  %Numero punti di controllo
Gamma = ones(m,n);

Xc = reshape(Xc,Nc,1,3);
N = reshape(N, Nc, 1, 3);

A = [];

%Butto via un po' di memoria per creare una ripetizione del primo versore

Ndummy = zeros(m,n,3);

for k = 1 : Nc

    Ndummy(:,:,1)=repmat(N(k,1,1),m,n);
    Ndummy(:,:,2)=repmat(N(k,1,2),m,n);
    Ndummy(:,:,3)=repmat(N(k,1,3),m,n);


    V = Velocita(Xa, Xc(k,1,:), Gamma);
  
    riga = dot(V,Ndummy,3);
    riga = reshape(riga, 1, Nc);
    A = [A; riga];
  

end  


