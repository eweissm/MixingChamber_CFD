function v = bcGhost_v(v, t)
    global yf;
    % xf: x-cordinates of cell faces
    
     %get size of v
    M = size(v,1)-2;
    N = size(v,2)-1; 
    
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

end