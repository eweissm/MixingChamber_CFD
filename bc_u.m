function u = bc_u(u, t)
    global xf yc;
    % xf: x-cordinates of cell faces
    % yc: y-cordinates of cell centers incl. ghost cells
    
     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir

    %left BC -- x= 0, y=(0,2) --> defining u(1.5, j) as we are node based
    for j = 1:N+2
        if yc(j) <= 0.25
            u(1,j) = 0;
        elseif ( (yc(j) > 0.25) && (yc(j) < 0.75) )
            u(1,j) = (-1/3)*(u(3,j) - 4*u(2,j)); %nuemann condition w/ 2nd order Forward diff
        else % if y >= .75
            u(1,j) = 0;
        end
    end
    
    %Right Bc -- x = 3, y=(0,2) --> defining u(M+3/2, j) as we are node based
    for j = 1:N+2
        z1 = yc(j)-.25;
        g1 = -(3/4)*6*z1*(1-z1);

        if yc(j) <= 0.25
            u(M+1,j) = 0;
        elseif ( (yc(j) > 0.25) && (yc(j) < 1.25) )
            u(M+1,j) = g1;
        else % if y >= 1.25
            u(M+1,j) = 0;
        end
    end
    
    %Bottom BC -- y= 0, x=(0,3)--> defining u(i, 1) as we are cell centered
    for i = 1:M+1
        z2 = 2*(xf(i)-1.5);
        g2 = 6*z2*(1-z2);

        if xf(i) <= 1.5
            u(i,1) = -u(i,2);
        elseif ( (xf(i) > 1.5) && (xf(i) < 2) )
            u(i,1) = 2*g2-u(i,2);
        else % if x >= 2
            u(i,1) = -u(i,2);
        end
    end
    
    %Top BC -- y= 2, x=(0,3)--> defining u(i, N+2) as we are cell centered
    for i = 1:M+1
        z3 = 2*(xf(i)-2);
        g3 = sqrt(3)*(3/2)*6*z3*(1-z3);

        if xf(i) <= 2
            u(i,N+2) = -u(i,N+1);
        elseif ( (xf(i) > 2) && (xf(i) < 2.5) )
            u(i,N+2) = 2*g3-u(i,N+1);
        else % if x >= 2.5
            u(i,N+2) = -u(i,N+1);
        end
    end

end