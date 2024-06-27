function Y = parabolic_CN2_2D_Y3(Y, QY, dt)
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

    %applying BCs to v at t_n+1/2
    Y = bc_Y3(Y, t+dt/2);
    
    %calculate d
    for j = 1:N+6
        for i = 2:M+5
            d(i,j) = dh*Y(i+1,j)+(1-2*dh)*Y(i,j)+dh*Y(i-1,j) + QY(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+dt
    [a, b, c, d] = bcCN2_Y3(a, b, c, d, t+(dt));

    %loop over horizontal slices solving the tridiagonal system
    for i = 2:M+5
        Y(i,2:N+5) = mySolveTriDiag(a(i,2:N+5), b(i,2:N+5), c(i,2:N+5), d(i,2:N+5));
    end

    %applying BCs to v at t_n+dt
    Y = bc_Y3(Y, t+(dt));
end