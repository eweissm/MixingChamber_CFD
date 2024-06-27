function v = parabolic_CN1_2D_v(v, Qv, dt)
    global Re t h;

    %get size of u
    M = size(v,1)-2; %node based in x dir
    N = size(v,2)-1; %cell centered in y dir
    
    %declare a, b, c and d
    a = ones(M+2, N+1);
    b = ones(M+2, N+1);
    c = ones(M+2, N+1);
    d = zeros(M+2, N+1);

    %calculate d1/ d2 (equidistant mesh in both x and y, so both are equal
    %--> calling it dh arbitartily
    dh = dt/(2*Re*h^2);
    
    %setting a, b and c
    a= -dh*a;
    b = (1+2*dh)*b;
    c = -dh*c;

    %applying BCs to v at t_n
    v = bc_v(v, t);
    
    %calculate d
    for j = 2:N
        for i = 1:M+2
            d(i,j) = dh*v(i,j+1)+(1-2*dh)*v(i,j)+dh*v(i,j-1) + Qv(i,j)*(dt/2);
        end
    end

    %apply implicit BCs for t_n+1/2
    [a, b, c, d] = bcCN1_v(a, b, c, d, t+(dt/2));

    %loop over horizontal slices solving the tridiagonal system
    for j = 2:N
        v(2:M+1,j) = mySolveTriDiag(a(2:M+1,j), b(2:M+1,j), c(2:M+1,j), d(2:M+1,j));
    end

    %applying BCs to v at t_n+1/2
    v = bc_v(v, t+(dt/2));
end