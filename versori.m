function [N,T] = versori(wing)
% estrae le normali al pannello 
%
    for i = 1:size(wing,1)-1
        for j = 1:size(wing,2)-1

            v13 = reshape(wing(i+1,j+1,:)-wing(i,j,:),3,1);
            v24 = reshape(wing(i+1,j,:)-wing(i,j+1,:),3,1);
            N(i,j,:) = cross(v24,v13)/norm(cross(v24,v13));
            T(i,j,:)=(v13+v24)/norm(v13+v24);
        end
    end
end