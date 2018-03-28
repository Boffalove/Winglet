function ml = naca_mean_line(NACA4, n_nodes, chord)
% genera la linea media di un NACA 4 cifre
%
    M = floor(NACA4/1000);
    P = floor((NACA4 - 1000*M)/100);
    XX = floor(NACA4 - 1000*M - 100*P);
    
    p = P/10;
    m = M/100;
    

%     x=chord*0.5*(1-cos(pi*linspace(0,1,n_nodes)));
    x = linspace(0,chord,n_nodes);
    x_lead = x(x<p*chord);
    x_trail = x(x>=p*chord);
    ml(1,:) = x;
    ml(2,:) = [m/p^2 * ( 2*p*x_lead/chord-(x_lead/chord).^2),...
               m/(1-p)^2 *(1-2*p + 2*p*x_trail/chord-(x_trail/chord).^2)]*chord;
end