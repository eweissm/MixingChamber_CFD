function kx = calckx(u)
    global h;

     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
    %perform 2d composite midpoint integration
    kx = 0;
    for j = 2:N+1
        for i = 1:M
            %calculate kx
            kx = kx + (1/2)*((u(i+1,j) + u(i, j) )/2)^2*h^2;
        end
   end 
end