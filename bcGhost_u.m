function u = bcGhost_u(u, t)
    global xf;
    % xf: x-cordinates of cell faces
    
     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
 %Bottom BC -- y= 0, x=(0,3)--> defining u(i, 1) as we are cell centered
    for i = 1:M+1
        z2 = 2*(xf(i)-1.5);
        g2 = 6*z2*(1-z2);

        if ( (xf(i) > 1.5) && (xf(i) < 2) )
            u(i,1) = 2*g2-u(i,2);
        else % if x >= 2
            u(i,1) = -u(i,2);
        end
    end
    
    %Top BC -- y= 2, x=(0,3)--> defining u(i, N+2) as we are cell centered
    for i = 1:M+1
        z3 = 2*(xf(i)-2);
        g3 = sqrt(3)*(3/2)*6*z3*(1-z3);

        if ( (xf(i) > 2) && (xf(i) < 2.5) )
            u(i,N+2) = 2*g3-u(i,N+1);
        else % if x >= 2.5
            u(i,N+2) = -u(i,N+1);
        end
    end

end