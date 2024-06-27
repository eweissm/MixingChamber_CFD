function [a, b, c, d] = bcCN2_Y3(a, b, c, d, t)

    global xc3;
    
    %get size of y
    M = size(a,1)-6; 
    N = size(a,2)-6; 

    %for each col
    for i = 1:M+6
        %apply bottom BC
        if((xc3(i) >= 1.5) && (xc3(i) <= 2) ) %inlet condition
           b(i,2:4) = b(i,2:4)-a(i,2:4);
           a(i,2:4) = 0;
        else %wall condition
           b(i,2:4) = b(i,2:4)+a(i,2:4);
           a(i,2:4) = 0;        
        end
        
        %apply top BC
        if((xc3 (i) >= 2) && (xc3(i) <= 2.5) ) %inlet condition  
           b(i,N+3:N+5) = b(i,N+3:N+5)-c(i,N+3:N+5);
           c(i,N+3:N+5) = 0;
        else %wall condition
           b(i,N+3:N+5) = b(i,N+3:N+5)+c(i,N+3:N+5);
           c(i,N+3:N+5) = 0;
        end
    end
end