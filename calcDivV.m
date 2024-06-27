function divV = calcDivV(u, v)
    global h;

     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
    %initialized divV withj zeros for all cells (including ghost cells)
    divV = zeros(M+2, N+2);

    for j = 2:N+1
        for i = 2:M+1
            %calculate 2nd order finite diff for divV
            divV(i,j) = ((u(i,j) - u(i-1,j) ) / (h)) +((v(i,j) - v(i,j-1) ) / (h)) ;

        end
    end
end