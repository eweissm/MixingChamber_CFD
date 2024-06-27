function du = calcdu(u)
    global h;
     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
    %initialized du
    du = zeros(M+2, N+2);

    for j = 2:N+1
        for i = 2:M+1
            %calculate 2nd order finite diff
            du(i,j) = (u(i,j) - u(i-1,j) ) / (h);

        end
    end
end