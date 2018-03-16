%% ala
clear all
close all

naca = 4412;
wing_root = 1;
wing_tip = 1;
wing_span = 10;
n_chord = 30;
n_span = 10;
wing_sweep = deg2rad(0);
wing_twist = deg2rad(0);

% winglet
w_let_root = wing_tip;
w_let_tip = .1;
height = 2;
Radius = .5;
n_height = 20;
cant = deg2rad(15);
sweep = deg2rad(15);
toe_out = deg2rad(-20);
up = 1; % 1 --> winglet verso l'alto
         % -1 --> winglet verso il basso
         
V_inf = zeros(n_chord-1,n_span-1,3);
V_inf(:,:,1) = -1;
V_inf(:,:,3) = .1;

dt = 1;
         
%% code
Wing = build_wing(wing_root,wing_tip,wing_span,n_chord,n_span,naca,...
                                                    wing_sweep,wing_twist);

% W_let = build_winglet(w_let_root,w_let_tip,height,Radius,n_chord,...
%                     n_height,cant,sweep,toe_out,naca,Wing,wing_twist,up);
% V_inf = zeros(n_chord-1,n_span-1+n_height-1,3);
% V_inf(:,:,1) = -1;
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
[~,A] = Induzione(Wing,p_controllo,N,ones(size(N,1),size(N,2)));
B = -dot(V_inf,N,3)';
b = B(:);

% reshape riempie per colonne quindi devo invertire gli indici e trasporre
gammad = A\b;
gamma = reshape(gammad,size(N,2),size(N,1))';

%% creo prima fila di pannelli di scia
V_wake = zeros(1,size(Wing,2),3);
for i = 1: size(Wing,2)
    [~,V_wake(1,i,:)] = Velocita(Wing,Wing(end,i,:),gamma); 
end
V_wake = repmat(-V_inf(1,1,:),1,size(Wing,2),1)+V_wake;
wake = [Wing(end,:,:);Wing(end,:,:)+V_wake*dt];
g_wake = gamma(end,:);

% figure(350)
% surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
% hold on
% axis equal
% plot3(wake(:,:,1),wake(:,:,2),wake(:,:,3),'r*')



%% ciclo temporale
for t = 1:50
    % risolvo sitema lineare
    V_sa = zeros(size(p_controllo));
    for i = 1: size(p_controllo,1)
        parfor j = 1:size(p_controllo,2)
            [~,V_sa(i,j,:)] = Velocita(wake,p_controllo(i,j,:),g_wake); 
        end
    end
    B = -dot(V_inf+V_sa,N,3)';
    b = B(:);
    gammad = A\b;
    gamma = reshape(gammad,size(N,2),size(N,1))';

    % muovo scia
    V_inf_wake = repmat(V_inf(1,1,:),size(wake,1),size(wake,2),1);
    
    wake = update_scia(Wing,gamma,wake,g_wake,dt)-V_inf_wake*dt;
    
    wake = [Wing(end,:,:); wake];
    g_wake = [gamma(end,:); g_wake];
    
end
