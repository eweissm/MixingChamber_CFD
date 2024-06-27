function [a, b, c, d] = bcCN1_Y3(a, b, c, d, t)

    global yc3;

    %get size of y
    M = size(a,1)-6; 
    N = size(a,2)-6; 

    %for each row
    for j = 1:N+6
        %apply left BC
        b(2:4,j) = b(2:4,j)+a(2:4,j);
        a(2:4,j) = 0;
        %c and d unchanged
        
        %apply right BC
        if ( (yc3(j) >= .25) && (yc3(j) <= 1.25) ) %inlet condition
            b(M+3:M+5, j) = b(M+3:M+5, j) -  c(M+3:M+5, j);
            d(M+3:M+5, j) = d(M+3:M+5, j) -  2* c(M+3:M+5, j);
            c(M+3:M+5, j) = 0;
            %a unchanged

        else %wall condition
            b(M+3:M+5, j) = b(M+3:M+5, j) +  c(M+3:M+5, j);
            c(M+3:M+5, j) = 0;
            %a and d unchanged
        end
       
    end
end