function f = calc_bdYdy_WENO_2D(Y,b)
    global h;
    
    %get size of Y
    M = size(Y,1)-6;
    N = size(Y,2)-6;

    %declare f --  Set Ghost cells to zero (this will not be changed)
    f = zeros(M+6,N+6);
    dp = zeros(M+6,N+6);
    dmdp = zeros(M+6,N+6);

    % get dp and dmdp
    dp(4:M+3,1:N+5) = Y(4:M+3,2:N+6) - Y(4:M+3, 1:N+5);
    dmdp(4:M+3,2:N+5) = Y(4:M+3,3:N+6) - 2*Y(4:M+3, 2:N+5) + Y(4:M+3, 1:N+4);

    %calculate f 
    for i = 4:M+3
        for j = 4:N+3

            %get WENO value --> will be added to dphidx
            if(b(i,j)>=0)
                WENO_val = psiWENO(dmdp(i, j-2)/h, dmdp(i, j-1)/h,dmdp(i, j)/h,dmdp(i, j+1)/h);
            else
                WENO_val = -psiWENO(dmdp(i, j+2)/h, dmdp(i, j+1)/h,dmdp(i, j)/h,dmdp(i, j-1)/h);
            end
            
            dphidx = (1/(12*h)) *(-dp(i, j-2)+ 7*dp(i, j-1) +7*dp(i, j) - dp(i, j+1)) - WENO_val;
            
            f(i,j) = b(i,j)*dphidx;
        
        end
    end

end