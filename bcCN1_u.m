function [a, b, c, d] = bcCN1_u(a, b, c, d, t)

    global yc;

    %get size of u
    M = size(a,1)-1; %node based in x dir
    N = size(a,2)-2; %cell centered in y dir

    %for each row
    for j = 2:N+1
        %apply left BC
        if((yc(j) > 0.25) && (yc(j) < 0.75) ) %outlet condition
            b(2, j) = b(2, j) + (4/3)*a(2, j);
            c(2, j) = c(2, j) - (1/3)*a(2, j);
            a(2, j) = 0;
            %d2 unchanged
        else %wall condition
            a(2, j) = 0;
            %d2, b2, c2 unchanged
        end
        
        %apply right BC
        if((yc(j) > 0.25) && (yc(j) < 1.25) ) %inlet condition
            z1 = yc(j)-.25;
            g1 = -(3/4)*6*z1*(1-z1);
            
            d(M, j) = d(M, j) - c(M, j)*g1;
            c(M, j) = 0;
            
            % b_M, a_M unchanged
        else %wall condition
            c(M, j) = 0;
            %d_M, b_M, a_M unchanged
        end
    end
end