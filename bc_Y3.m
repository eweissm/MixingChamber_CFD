function Y = bc_Y3(Y,t)
    global xc3 yc3;

    %get size of y
    M = size(Y,1)-6;
    N = size(Y,2)-6;

    for j = 1:N+6 %run through all columns

        %left BC
        Y(3,j) = Y(4,j);
        Y(2,j) = Y(5,j);
        Y(1,j) = Y(6,j);

        %Right Bc 
        if ( (yc3(j) >= .25) && (yc3(j) <= 1.25) ) %inlet condition
            Y(M+4,j) = 2-Y(M+3,j);
            Y(M+5,j) = 2-Y(M+2,j);
            Y(M+6,j) = 2-Y(M+1,j);
        else                                         %wall condition
            Y(M+4,j) = Y(M+3,j);
            Y(M+5,j) = Y(M+2,j);
            Y(M+6,j) = Y(M+1,j);
        end

    end
    
    for i = 1:M+6

        %Bottom BC 
       if ( (xc3(i) >= 1.5) && (xc3(i) <= 2) ) %inlet condition
            Y(i,3) = -Y(i,4);
            Y(i,2) = -Y(i,5);
            Y(i,1) = -Y(i,6);
       else                                %wall condition
            Y(i,3) = Y(i,4);
            Y(i,2) = Y(i,5);
            Y(i,1) = Y(i,6);
       end

       
       %Top BC 
       if ( (xc3(i) >= 2) && (xc3(i) <= 2.5) ) %inlet condition
           Y(i,N+4) = -Y(i,N+3);
           Y(i,N+5) = -Y(i,N+2);
           Y(i,N+6) = -Y(i,N+1);
       else                                  %wall condition
           Y(i,N+4) = Y(i,N+3);
           Y(i,N+5) = Y(i,N+2);
           Y(i,N+6) = Y(i,N+1);
       end

    end

end