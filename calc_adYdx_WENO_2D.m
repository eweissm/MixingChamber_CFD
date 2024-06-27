function f = calc_adYdx_WENO_2D(Y,a)
    global h;
    
    %get size of Y
    M = size(Y,1)-6;
    N = size(Y,2)-6;

    %declare f --  Set Ghost cells to zero (this will not be changed)
    f = zeros(M+6,N+6);
    dp = zeros(M+6,N+6);
    dmdp = zeros(M+6,N+6);
    
    % get dp and dmdp
    dp(1:M+5,4:N+3) = Y(2:M+6,4:N+3) - Y(1:M+5,4:N+3);
    dmdp(2:M+5, 4:N+3) = Y(3:M+6,4:N+3) - 2*Y(2:M+5,4:N+3) + Y(1:M+4,4:N+3);
    
    %calculate f
    for j = 4:N+3
        for i = 4:M+3

            %get WENO value --> will be added to dphidx
            if(a(i,j)>=0)
                WENO_val = psiWENO(dmdp(i-2,j)/h, dmdp(i-1,j)/h,dmdp(i,j)/h,dmdp(i+1,j)/h);
            else
                WENO_val = -psiWENO(dmdp(i+2,j)/h, dmdp(i+1,j)/h,dmdp(i,j)/h,dmdp(i-1,j)/h);
            end
            
            dphidx = (1/(12*h)) *(-dp(i-2,j)+ 7*dp(i-1,j) +7*dp(i,j) - dp(i+1,j)) - WENO_val;
            
            f(i,j) = a(i,j)*dphidx;
        
        end
    end

end