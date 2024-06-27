function [S] = calcS3(Y)
    global Lx Ly h;


         %get size of Mesh
    M = size(Y,1)-6; 
    N = size(Y,2)-6;
    
    S = 0;

    for j = 4:N+3
        for i = 4:M+3
            %calculate S
            S = S + Y(i,j)*(1-Y(i,j));
        end
    end 
    

    S = S*(h^2) / (Lx*Ly);
end