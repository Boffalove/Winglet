%% ala
clear all
close all

naca = 4412;
wing_root = 1;
wing_tip = 1;
wing_span = 10;
n_chord = 10;
n_span = 10;
wing_sweep = deg2rad(0);
wing_twist = deg2rad(0);

% winglet
w_let_root = wing_tip;
w_let_tip = .1;
height = 2;
Radius = .25;
n_height = 20;
cant = deg2rad(15);
sweep = deg2rad(15);
toe_out = deg2rad(-20);
up = 1; % 1 --> winglet verso l'alto
         % -1 --> winglet verso il basso
         
V_inf = zeros(n_chord-1,n_span-1,3);
V_inf(:,:,1) = -1;
V_inf(:,:,3) = 0.1;
         
%% code
Wing = build_wing(wing_root,wing_tip,wing_span,n_chord,n_span,naca,...
                                                    wing_sweep,wing_twist);

% W_let = build_winglet(w_let_root,w_let_tip,height,Radius,n_chord,...
%                     n_height,cant,sweep,toe_out,naca,Wing,wing_twist,up);
% V_inf = zeros(n_chord-1,n_span-1+n_height-1,3);
% V_inf(:,:,1) = 1;
% V_inf(:,:,3) = .1;
% Wing = assemble_wing(Wing,W_let);                                                

[N,T] = versori(Wing);
p_controllo = collocazione(Wing);

%% plot ala
% 
% figure(100)
% surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
% xlabel('chord')
% ylabel('span')
% axis equal
% hold on
% plot3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),'r*')
% quiver3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),...
%                                           N(:,:,1),N(:,:,2),N(:,:,3),'k')
% quiver3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),...
%                                           T(:,:,1),T(:,:,2),T(:,:,3),'y')

%% risoluzione sistema lineare
tic
[~,A] = Induzione(Wing,p_controllo,N,ones(size(N,1),size(N,2)));
toc
B = -dot(V_inf,N,3)';
b = B(:);

% reshape riempie per colonne quindi devo invertire gli indici e trasporre
gammad = A\b;
gamma = reshape(gammad,size(N,2),size(N,1))';

%% plot gamma
figure(200)
surf(gamma)
% superficie pannelli sempre uguale solo se profilo simmetrico e spaziatura 
% uniforme(non coseni);
S_pann = sqrt(sum((Wing(1,2,:)-Wing(1,1,:)).^2,3))*...
                                sqrt(sum((Wing(2,1,:)-Wing(1,1,:)).^2,3));
                            
G_tot = sum(sum(gamma))*S_pann;
% pressione
V = Induzione(Wing,p_controllo,N,gamma);
cp = 1-(sqrt(sum((V+V_inf).^2,3)).^2./sqrt(sum(V_inf.^2,3))).^2;
figure(300)
surf(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),cp)
axis equal
colorbar


%%
[lines,g_lines] = build_lines(Wing,gamma);
for i = 1:size(p_controllo,1)
    for j = 1: size(p_controllo,2)
        Vl(i,j,:) = sum(induced_line_speed(lines,g_lines,p_controllo(i,j,:)),1);
    end
end
cpl = 1-(sqrt(sum((Vl+V_inf).^2,3)).^2./sqrt(sum(V_inf.^2,3))).^2;
figure(400)
surf(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),cpl)
axis equal
colorbar