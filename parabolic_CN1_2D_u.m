function u = parabolic_CN1_2D_u(u, Qu, dt)
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

    %applying BCs to u at t_n
    u = bc_u(u, t);
    
    %calculate d
    for j = 2:N+1
        for i = 1:M+1
            d(i,j) = dh*u(i,j+1)+(1-2*dh)*u(i,j)+dh*u(i,j-1) + Qu(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+1/2
    [a, b, c, d] = bcCN1_u(a, b, c, d, t+(dt/2));

    %loop over horizontal slices solving the tridiagonal system
    for j = 2:N+1
        u(2:M,j) = mySolveTriDiag(a(2:M,j), b(2:M,j), c(2:M,j), d(2:M,j));
    end

    %applying BCs to u at t_n+1/2
    u = bc_u(u, t+(dt/2));
end