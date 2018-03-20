function [u_wake] = update_scia(wing,g_wing,wake,g_wake,dt)
% aggiorna la scia di un vortex lattice instazionario 
% 
    V = zeros(size(wake,1),size(wake,2),3);
    
    for i = 1:size(wake,1)
        parfor j = 1:size(wake,2)
            % ala su scia
            [~ , vv_as]=Velocita(wing,wake(i,j,:),g_wing);
            % scia su scia 
            [~,vv_ss]=Velocita(wake,wake(i,j,:),g_wake);
            V(i,j,:) = vv_as+vv_ss; 
        end
    end
    u_wake = wake+V*dt;
end

