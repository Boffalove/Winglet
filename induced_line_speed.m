function V = induced_line_speed(line,gamma,point)
% Line = [ fine_vortice  inizio_vortice]
% calcola la velocit� indotta da una linea vorticosa avente intensit� gamma
% nel punto / nei punti descritti da point
%
% line 1,2,3 o matrice colonna
% gamma = 1 o vettore colonna
% point 1,1,3 o n,n,3

    r1 = point-line(:,2,:);
    r2 = point-line(:,1,:);

    % if sqrt(sum(r1.^2,3)) == 0 || sqrt(sum(r2.^2,3)) == 0
    % V = zeros(1,1,3);
    % else
    [riga1,colonna1] = find(sqrt(sum(r1.^2,3))==0);
    [riga2,colonna2] = find(sqrt(sum(r2.^2,3))==0);

    num =(sqrt(sum(r1.^2,3))+sqrt(sum(r2.^2,3))).*cross(r1,r2,3);
    den = sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3)).*...
          (sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3))+dot(r1,r2,3));
    % correzione per tenere conto del core viscoso
    k = 0;1e-3; %parametro modificabile
    corr1 = k^2*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)).^2)./ ...
                sqrt(sum(cross(r1,r2,3).^2,3)).^2;
    corr2 = 2*k*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)))./ ...
                sqrt(sum(cross(r1,r2,3).^2,3));
    corr = 1 + corr1 + corr2;
    
    V = gamma./(4*pi)*num./(den.*corr);
    
    if ~isempty(riga1)
        V(riga1,colonna1,:) = zeros(size(riga1,2),size(colonna1,2),3);
    end
    if ~isempty(riga2)
        V(riga2,colonna2,:) = zeros(size(riga2,2),size(colonna2,2),3);
    end   
end