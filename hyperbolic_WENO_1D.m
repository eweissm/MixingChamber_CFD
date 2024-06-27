function Y = hyperbolic_WENO_1D(Y, a0, a1, a2, dt)
    global t;
    
    %get size of Y
    M = size(Y,1)-6;

    Y1 = zeros(M+6,1);
    Y2 = zeros(M+6,1);

    %assuming BC have previously been applied to Y

    %step 1
    dphi_n_dx = calc_adYdx_WENO(Y, a0);
    Y1 = Y - dt * dphi_n_dx;
    Y1 = bc3(Y1, t+dt);
    
    %step 2
    dphi_1_dx = calc_adYdx_WENO(Y1, a1);
    Y2 = Y1 +(3/4)*dt*dphi_n_dx - (1/4)*dt*dphi_1_dx;
    Y2 = bc3(Y2, t+(1/2)*dt);

    %step 3
    Y = Y2+(1/12)*dt*dphi_n_dx + (1/12)*dt*dphi_1_dx - (2/3)*dt*calc_adYdx_WENO(Y2, a2);
    Y = bc3(Y, t+dt);

end
