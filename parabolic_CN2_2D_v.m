function v = parabolic_CN2_2D_v(v, Qv, dt)
    global Re t h;

    %get size of v
    M = size(v,1)-2; % x dir
    N = size(v,2)-1; % y dir
    
    %declare a, b, c and d
    a = ones(M+2, N+1);
    b = ones(M+2, N+1);
    c = ones(M+2, N+1);
    d = zeros(M+2, N+1);

    %calculate d1/d2 (equidistant mesh in both x and y, so both are equal
    %--> calling it dh arbitartily
    dh = dt/(2*Re*h^2);
    
    %setting a, b and c
    a= -dh*a;
    b = (1+2*dh)*b;
    c = -dh*c;

    %applying BCs to v at t_n (t_n should be t+dt)
    v = bc_v(v, t+dt/2);
    
    %calculate d
    for j = 1:N+1
        for i = 2:M+1
            d(i,j) = dh*v(i+1,j)+(1-2*dh)*v(i,j)+dh*v(i-1,j) + Qv(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+1/2
    [a, b, c, d] = bcCN2_v(a, b, c, d, t+dt);

    %loop over horizontal slices solving the tridiagonal system
    for i = 2:M+1
        v(i,2:N) = mySolveTriDiag(a(i,2:N), b(i,2:N), c(i,2:N), d(i,2:N));
    end

    %applying BCs to v at t_n+1/2
    v = bc_v(v, t+(dt));
end