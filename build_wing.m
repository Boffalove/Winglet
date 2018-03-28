function wing = build_wing(root_chord, tip_chord, span, n_chord,...
                           n_span, NACAs, sweep, twist, varargin)
% costruisce 3 matrici contenenti la posizione dei nodi di discretizzazione
% dell'ala
    if length(NACAs) == 1
        NACAs = repmat(NACAs,1,n_span);
    end
    chord_vect = linspace(root_chord, tip_chord, n_span);
%     yy=span*0.5*(1-cos(pi*linspace(0,1,n_span)));
    
    yy = linspace(0,span, n_span);
    
    for i = 1:n_span
        ml = naca_mean_line(NACAs(i), n_chord, chord_vect(i));
        wing(:,i,1) = ml(1,:);
        wing(:,i,3) = ml(2,:);
    end
    
    % aggiungo il twist
    if nargin >= 9
        old_twist = varargin{2};
        twist_vect = linspace(old_twist,old_twist+twist,n_span);
    else
        twist_vect = linspace(0,twist,n_span);
    end
    wing(:,:,1) = wing(:,:,1).*repmat(cos(twist_vect),n_chord,1) +...
                  wing(:,:,3).*repmat(sin(twist_vect),n_chord,1);
    wing(:,:,3) = -wing(:,:,1).*repmat(sin(twist_vect),n_chord,1) +...
                  wing(:,:,3).*repmat(cos(twist_vect),n_chord,1);
    % allineo il punto medio del profilo 
    wing(:,:,1) = wing(:,:,1) + repmat((root_chord -chord_vect)/2,n_chord,1);
    % aggiungo lo sweep 
    wing(:,:,1) = wing(:,:,1) + yy*sin(sweep);
    wing(:,:,2) = repmat(yy,n_chord,1).*cos(sweep);

    if nargin >= 9
        old_wing = varargin{1};
        wing(:,:,1) = wing(:,:,1)+old_wing(1,end,1);
        wing(:,:,2) = wing(:,:,2)+old_wing(1,end,2);

%         wing(:,:,3) = wing(:,:,3);

    end
end