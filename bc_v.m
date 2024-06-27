function v = bc_v(v, t)
    global xc yf;
    % xc: x-cordinates of cell centers incl. ghost cells
    % yf: y-cordinates of cell faces
    
    %get size of v
    M = size(v,1)-2; % in x dir
    N = size(v,2)-1;  %y dir

    
    for j = 1:N+1

        %left BC
        if ( (yf(j) >= 0.25) && (yf(j) <= 0.75) ) %outlet condition
            v(1,j) = v(2,j);
        else                                    %wall condition
            v(1,j) = -v(2,j); 
        end

        %Right Bc 
        v(M+2,j) = -v(M+1,j);
    end
    
    for i = 1:M+2

        %Bottom BC 
       if ( (xc(i) >= 1.5) && (xc(i) <= 2) ) %inlet condition
            z2 = 2*(xc(i)-1.5);
            g2 = sqrt(3)*6*z2*(1-z2);

            v(i,1) = g2;
        else                                %wall condition
            v(i,1) = 0;
       end

       
       %Top BC 
       if ( (xc(i) >= 2) && (xc(i) <= 2.5) ) %inlet condition
            z3 = 2*(xc(i)-2);
            g3 = -9*z3*(1-z3);
            v(i,N+1) = g3;
       else                                  %wall condition
            v(i,N+1) = 0;
       end

    end

end