%Prova funzioni-------------_>tempi

clear; close all; clc;

M = 10;               % Numero pannelli in corda
N = 10 : 5: 200;      % Numero pannelli in apertura

tempi = [];


for n = N
  
  Xa = rand(M+1, n+1, 3);
  Xc = rand(M, n, 3);
  Normali = rand(M, n, 3);
  
  tic
  A = CalcolaMatrice(Xa, Xc, Normali);
  t = toc;
  
  tempi = [tempi; t];
  
end

figure(1)
plot(N + M*ones(1, length(N)), tempi, 'o-','linewidth',2);
grid on
xlabel('Numero pannelli')
ylabel('Tempo calcolo [s]')

