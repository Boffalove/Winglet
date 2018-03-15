function [V,A] = Induzione(Xa,Xb,N,Gamma)

%Calcola la velocit√† indotta nei punti della griglia Xb
%dalla griglia di vortici Xa. La griglia di vortici √® 
%formata da anelli vorticosi: ogni elemento (3D) di Xa
%corrisponde ad un vertice dell'anello. La circolazione
%di ogni anello √® contenuta nella matrice Gamma.
%La funzione calcola anche la matrice A:
%
%     A(i,j) = velocit√†_anello_j * N_i
%
%     dove N_i √® il versore normale al pannello i-esimo 
%
%La matrice (3D) N contiene i versori normali dei pannelli 
%i cui punti di controllo sono gli elementi della griglia Xb

% Xa : (ma,na,3)
% Xb : (mb,nb,3)
% N  : (mb,nb,3)
% Gamma : (ma,na)
% V :  (mb,nb,3)
% A :  (mb,nb)

%!!!!!!!E' talmente lenta da essere imbarazzante!!!!!!!!!!!
%       CRIBBIOOOOOOO!!!!!!


ma = size(Xa,1);
na = size(Xa,2);
mb = size(Xb,1);
nb = size(Xb,2);

ring = zeros(2,2,3);
V = zeros(mb,nb,3);
A = [];
for i = 1 : ma - 1
   for j = 1 : na - 1
      
%       ring = [Xa(i,j,:), Xa(i,j+1,:);
%               Xa(i+1,j,:), Xa(i+1,j+1,:)];
          
      ring = Xa(i:i+1,j:j+1,:);

%       V = V + induced_ring_speed(ring, Gamma(i,j), Xb);
      
      Vij =  induced_ring_speed(ring, Gamma(i,j), Xb); % velocit‡ indotta
                              % da un vortice su tutti i punti di controllo
      V = V+Vij;
      
      Vijd = (dot(Vij,N,3))'; % mi serve sotto
      % devo riordinare perchË cosi ho sulla colonna di A la velocit‡ 
      % indotta dal vortice corrispondente alla colonna stessa
      A(:,end+1)=Vijd(:)';
      
      
   end
end

% A = sum(V.*N, 3);% A = dot(V,N,3);

end
