function u = parabolic_CN2_2D_u(u, Qu, dt)
    global Re t h;

    %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
    %declare a, b, c and d
    a = ones(M+1, N+2);
    b = ones(M+1, N+2);
    c = ones(M+1, N+2);
    d = zeros(M+1, N+2);

    %calculate d1/ d2 (equidistant mesh in both x and y, so both are equal
    %--> calling it dh arbitartily
    dh = dt/(2*Re*h^2);
    
    %setting a, b and c
    a= -dh*a;
    b = (1+2*dh)*b;
    c = -dh*c;

    %applying BCs to u at t_n+1/2
    u = bc_u(u, t+dt/2);
    
    %calculate d
    for j = 1:N+2
        for i = 2:M
            d(i,j) = dh*u(i+1,j)+(1-2*dh)*u(i,j)+dh*u(i-1,j) + Qu(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+dt
    [a, b, c, d] = bcCN2_u(a, b, c, d, t+(dt));

    %loop over horizontal slices solving the tridiagonal system
    for i = 2:M
        u(i,2:N+1) = mySolveTriDiag(a(i,2:N+1), b(i,2:N+1), c(i,2:N+1), d(i,2:N+1));
    end

    %applying BCs to u at t_n+dt
    u = bc_u(u, t+(dt));
end