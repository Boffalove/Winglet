function V = induced(vortex,gamma,points,N)
% calcola la velocità indotta da una griglia di vortici a staffa descritta
% dai punti contenuti in vortex(m,n,3) aventi intensità gamma(m-1,n-1)
% nei punti points(k,l,3) 
% se numel(gamma)==1 allora calcola i coefficienti d'induzione e allora
% [l,k] = [m-1,n-1]
% V(k,l,3) se numel gamma~=1
% V(m-1*n-1,m-1*n-1) se numel(gamma)==1

    % preallocazione memoria
    if numel(gamma)==1 % voglio i coefficienti d'induzione
        V = [];
    else % voglio la velocità indotta
        V = zeros(size(points));
    end
    
    % numero di punti in apertura
    pp = size(points,2);
    
    % cilco sull'apertura
    for j = 1:size(points,2)
        % ciclo sulla corda
        for i = 1:size(points,1)
            % velcotià indotta dagli spanwise
            vs = span_induced(vortex,gamma,points(i,j,:));
            % velocità indotta dai chordwise
            vc = chord_induced(vortex,gamma,points(i,j,:));

            if numel(gamma)==1 % voglio una matrice 2D delle velocità normali

                dummy = vs+vc;
                dummy = reshape(dummy,1,size(dummy,1)*size(dummy,2),3);
                V = [V;dot(dummy,repmat(N(i,j,:), ...
                                                   1,size(dummy,2),1),3)];
        
            else % volgio una matrice 3D di velocità indotte
                
                V(i,j,:) = vs+vc;
            end
            
        end
    end
end