function Y = parabolic_CN1_2D_Y3(Y, QY, dt)
    global Re Sc t h;

   %get size of y
    M = size(Y,1)-6; 
    N = size(Y,2)-6; 

    
    %declare a, b, c and d
    a = ones(M+6, N+6);
    b = ones(M+6, N+6);
    c = ones(M+6, N+6);
    d = zeros(M+6, N+6);

    %calculate d1/ d2 (equidistant mesh in both x and y, so both are equal
    %--> calling it dh arbitartily
    dh = dt/(2*Re*Sc*h^2);
    
    %setting a, b and c
    a= -dh*a;
    b = (1+2*dh)*b;
    c = -dh*c;

    %applying BCs to v at t_n
    Y = bc_Y3(Y, t);
    
    %calculate d
    for j = 2:N+5
        for i = 1:M+6
            d(i,j) = dh*Y(i,j+1)+(1-2*dh)*Y(i,j)+dh*Y(i,j-1) + QY(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+1/2
    [a, b, c, d] = bcCN1_Y3(a, b, c, d, t+(dt/2));

    %loop over horizontal slices solving the tridiagonal system
    for j = 1:N+6
        Y(2:M+5,j) = mySolveTriDiag(a(2:M+5,j), b(2:M+5,j), c(2:M+5,j), d(2:M+5,j));
    end

    %applying BCs to v at t_n+1/2
     Y = bc_Y3(Y, t+(dt/2));
end