function VV = chord_induced(vortex,gamma,point)
% calcola la velocità indotta dai segmenti allineati alla linea media del
% profilo dei vortici ad ANELLO descritti dalla matrice 3d di punti
% vortex(m,n,3) aventi intensità gamma nel punto point(1,1,3)
% tiene conto anche dei filamenti di 'scia' per questo serve la direzione
% lungo la quale sono allineati dir(1,1,3)
% se gamma ha un solo elemento allora VV è la matrice che descrive
% l'induzione di tutti i vortici in quel punto considerando intensità
% unitaria
% se gamma(m-1,n-1) allora VV è la velocità indotta da tutti i vortici nel
% punto scelto VV(1,1,3);

    % assegno il valore di gamma corrispondente al segmento
    if numel(gamma ==1) % coefficienti d'influenza
       gamma_c = ones(size(vortex,1)-1,size(vortex,2));
    else
        % per ogni segmento vorticoso assegno il valore di gamma
        % corrispondente
        % per ogni riga di vortex il segmento j-esimo avra intensità pari a
        %  gamma_i - gamma_(i+1) dove i è l'indice dei vortici ad anello
        gamma_c = gamma; gamma_c(:,2:end) = +gamma(:,2:end)-gamma(:,1:end-1);
        gamma_c(:,end+1) = gamma(:,end); % numel(segmenti_riga) = numel(vortici_riga)+1
    end
    
    % estraggo i segmenti, presi tutti nella stessa direzione (verso il
    % bordo d'attacco
    % punto finale
    line_e = reshape(vortex(1:end-1,:,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)),1,3);
    %punto iniziale
    line_s = reshape(vortex(2:end,:,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)),1,3);
                                    
    gamma_c = reshape(gamma_c,numel(gamma_c),1);
    
    % costriusco il segmento
    line = [line_e,line_s];
    
    
    % calcolo della velocità indotta dei segmenti
    
    r1 = point-line(:,2,:);
    r2 = point-line(:,1,:);


    %cerco h == 0 per imporre poi velocità nulla
    [riga1,colonna1] = find(sqrt(sum(r1.^2,3))==0);
    [riga2,colonna2] = find(sqrt(sum(r2.^2,3))==0);
    [riga3,colonna3] = find(sqrt(sum(cross(r1,r2,3).^2,3))==0);
    
    
    % numeratore e denominatore per Biot-Savart
    num =(sqrt(sum(r1.^2,3))+sqrt(sum(r2.^2,3))).*cross(r1,r2,3);
    den = sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3)).*...
          (sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3))+dot(r1,r2,3));
      
    % correzione per tenere conto del core viscoso
    % deriva dall'aggiunta a denominatore di un termine(raggio del core
    % viscoso) che serve per "simulare" l'effetto del core allontanado il
    % punto in modo artificiale
    
    k = 0; %parametro modificabile
    
    % sviluppato analiticamente per ottenere una forma di più facile
    % implementazione numerica
    if k~=0
        corr1 = k^2*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)).^2)./ ...
                    sqrt(sum(cross(r1,r2,3).^2,3)).^2;
        corr2 = 2*k*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)))./ ...
                    sqrt(sum(cross(r1,r2,3).^2,3));

        corr = 1 + corr1 + corr2;
    else 
        corr = 1;
    end
    
    V = gamma_c/(4*pi).*num./(den.*corr);   
    
    % assegno nulle le velocità dove h=0
    if ~isempty(riga1)
        V(riga1,colonna1,:) = zeros(size(riga1,1),size(colonna1,1),3);
    end
    if ~isempty(riga2)
        V(riga2,colonna2,:) = zeros(size(riga2,1),size(colonna2,1),3);
    end   
    if ~isempty(riga3)
        V(riga3,colonna3,:) = zeros(size(riga3,1),size(colonna3,1),3);
    end  
    
    if numel(gamma)==1 %costruisco la matrice dei coefficienti d'induzione
                       % 3D
       V = reshape(V,size(vortex,1)-1,size(vortex,2),3);
       VV = ricostruzione_ind_vortice(V);
      
    else % restituisco la velocità indotta dai segmenti chordwise
       V = reshape(V,size(vortex,1)-1,size(vortex,2),3);
       VV = ricostruzione_ind_vortice(V);
       VV = sum(sum(VV));
    end
    
end


%=======================================================================
% funzioni usate per ricostruire la velocità indotta dalla matrice di
% vortici ad anello
function VV = ricostruzione_ind_vortice(V)
        VV =  V(:,1:end-1,:) - V(:,2:end,:);
%     for i = 1:size(V,1)
%         VV(i,:,:) = V
%         % primo vortice ad anello, tutte le righe prima e seconda colonna
%         % è l'unico ad essere gia giusto qunidi sommo il filamento di 
%         % sinistra e sottraggo quello di destra
%         VV(i,1,:) = +sum(V(i:end,1,:),1) - sum(V(i:end,2,:),1);
%         % tutti gli altri sono presi al contrario quindi sommo il vortice
%         % di sinistra e sottraggo quello di destra
%         VV(i,2:size(V,2)-1,:) = -sum(V(i:end,3:end,:),1) + ...
%                                                  sum(V(i:end,2:end-1,:),1);
%     end
end