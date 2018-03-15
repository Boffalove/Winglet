function V = induced_ring_speed(ring,gamma,point)
% calcola la velocità indotta da un anello vorticoso avente intensità gamma
% nel punto / nei punti descritti da point
%
% ring=2,2,3 
% gamma = 1
% point 1,1,3 o n,n,3

%   line1 = [ring(1,2,:) ring(1,1,:)]; % orizzontale sopra
%   line2 = [ring(2,2,:) ring(1,2,:)]; % verticale destra
%   line3 = [ring(2,1,:) ring(2,2,:)]; % orizzontale sotto
%   line4 = [ring(1,1,:) ring(2,1,:)]; % verticale sinistra

  line_v = [ring(1,2,:) ring(1,1,:); % orizzontale sopra
          ring(2,2,:) ring(1,2,:); % verticale destra
          ring(2,1,:) ring(2,2,:); % orizzontale sotto
          ring(1,1,:) ring(2,1,:)]; % verticale sinistra
%   r1 = point-line1(1,1,:);
%   r2 = point-line1(1,2,:);
  for i = 1:4
      line = line_v(i,:,:);
%       r1 = point-line_v(i,2,:);
%       r2 = point-line_v(i,1,:);
%   
%       num =(sqrt(sum(r1.^2,3))+sqrt(sum(r2.^2,3))).*cross(r1,r2,3);
%       den = sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3)).*...
%                  (sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3))+dot(r1,r2,3));
      if i == 1
%           V = num./den;
          V = induced_line_speed(line,gamma,point);
      else
%           V = V + num./den;
          V = V + induced_line_speed(line,gamma,point);
      end
  end
  
end