function HY = hyperbolic_Y_WENO_2D(Y, u, v, u1, v1, dt)
    global t
    
    M = size(Y,1)-6;
    N = size(Y,2)-6;

    a0 = zeros(M+6,N+6);
    a1 = zeros(M+6,N+6);
    b0 = zeros(M+6,N+6);
    b1 = zeros(M+6,N+6);

    %Generate a0, a1, a2, b0, b1 and b2 at cell centers
   
    a0(4:M+3,4:N+3) = (u(1:M, 2:N+1) + u(2:M+1,2:N+1))/2;
    a1(4:M+3,4:N+3) = (u1(1:M, 2:N+1) + u1(2:M+1,2:N+1))/2;
    b0(4:M+3,4:N+3) = (v(2:M+1, 1:N) + v(2:M+1, 2:N+1))/2;
    b1(4:M+3,4:N+3)= (v1(2:M+1, 1:N) + v1(2:M+1, 2:N+1))/2;
   
    a2 =(a0+a1) / 2;
    b2 =(b0+b1) / 2;

    %assuming BC have previously been applied to Y
    %step 1
    dphi_n_dx = calc_adYdx_WENO_2D(Y, a0);
    dphi_n_dy = calc_bdYdy_WENO_2D(Y, b0);

    Y1 = Y - dt * (dphi_n_dx + dphi_n_dy);  
    Y1 = bc_Y3(Y1, t+dt); 
    
    %step 2
    dphi_1_dx = calc_adYdx_WENO_2D(Y1, a1);
    dphi_1_dy = calc_bdYdy_WENO_2D(Y1, b1);

    Y2 = Y1 +(3/4)*dt*(dphi_n_dx + dphi_n_dy) - (1/4)*dt*(dphi_1_dx + dphi_1_dy);
    Y2 = bc_Y3(Y2, t+(1/2)*dt);

    %step 3
    dphi_2_dx=calc_adYdx_WENO_2D(Y2, a2);
    dphi_2_dy=calc_bdYdy_WENO_2D(Y2, b2);

    Ys = Y2+(1/12)*dt*(dphi_n_dx + dphi_n_dy) + (1/12)*dt*(dphi_1_dx + dphi_1_dy) - (2/3)*dt*(dphi_2_dx + dphi_2_dy);
    Ys = bc_Y3(Ys, t+dt);

    %calculate HY
    HY = (Ys - Y)/dt;
    
end