function k = calck(u, v)
    global h;

     %get size of Mesh
    M = size(u,1)-1; 
    N = size(u,2)-2;
    
    %perform 2d composite midpoint integration
    k = 0;
    for j = 2:N+1
        for i = 1:M
            %calculate kx
            k = k + ((u(i+1,j) + u(i, j) )/2)^2;
        end
    end 

    for j = 1:N
        for i = 2:M+1
            %calculate kx
            k = k + ((v(i,j+1) + v(i, j) )/2)^2;
        end
    end 
    

    k = 0.5*k*h^2;
end