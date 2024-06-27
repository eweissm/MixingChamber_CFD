clear
clear all
clf

global h CFL Re Sc xf xc yf yc yc3 xc3 Lx Ly t;

%Declare some constants
M=192;
N=128;
CFL = 0.25;
Re = 20;
Sc = 2;
Lx=3;
Ly=2;
times_of_interest = 0:.25:10;
Contours_times_of_interest = 4:.25:6;

x_range = [0,Lx];
y_range = [0,Ly];

%h is equal in both x and y
h = (x_range(2)-x_range(1))/(M);

%poissons solver controls:
nIterMax = 100;
epsilon = (h)/100; %Q is first order in time, so epsilon needs only be less than h

%lets define xf, yf, yc, xc, sc3 and yc3;
xf = linspace(x_range(1),x_range(2),M+1)';
yf = linspace(y_range(1),y_range(2),N+1)';
xc = linspace(x_range(1)-h/2,x_range(1)+(x_range(2)-x_range(1))+h/2,M+2)';
yc = linspace(y_range(1)-h/2,y_range(1)+(y_range(2)-y_range(1))+h/2,N+2)';
yc3 = linspace(y_range(1)-5*h/2,y_range(1)+(y_range(2)-y_range(1))+5*h/2,N+6)';
xc3 = linspace(x_range(1)-5*h/2,x_range(1)+(x_range(2)-x_range(1))+5*h/2,M+6)';

%lets define some variables we will be storing
%we do not know the length of the vectore before hand. Thus we will
%preallocate memory in blocks and add blocks as needed
BlockSize = 1000;% size increments added to vector length
currentSize = BlockSize; %current length of vectors
k_values = zeros(BlockSize,1);
S_values = zeros(BlockSize,1);
Time_values = zeros(BlockSize,1);

%declare t
t = 0;

%applying initial conditons at t=0--------------------------
    % 1. declare u, v and Y  --> set IC to 0
    u = zeros(M+1, N+2);
    v = zeros(M+2, N+1);
    Y = zeros(M+6, N+6);
    
    % 2. apply BCs to u, v and Y at T=0
    u = bc_u(u, t);
    v = bc_v(v, t);
    Y = bc_Y3(Y, t);
    
    %Adams-bashforth requires hyperbolic terms from previous step -->
    %initialize these previous terms for n-1 using values at n
    dt  = calcDtBurgers561(t, times_of_interest(2), u, v);
    HY_prev = hyperbolic_Y_WENO_2D(Y, u, v, u, v, dt); 
    [Hu_prev,Hv_prev] = hyperbolic_uv_2D(u,v);
    
    % 3. correct outlet velocities to ensure volume conservation
    %[u, v] = correctOutlet(u,v);
    
    % 4. calculate right hand side of Poisson equation using ∆t = 1
    % dt = 1;
    % f = (1/dt)*calcDivV(u, v);
    
    % 5. solve Poisson equation for Lagrange multiplier using Neumann boundary conditions
    % phi = ones(M+2, N+2); %inialized phi = 1
    % phi = bcGS(phi);
    % [phi, Linf, iter] = myPoisson(phi, f, h, nIterMax, epsilon);
    
    % 6. project velocities using Lagrange multiplier using ∆t = 1
    % [u, v]= projectV(u, v, phi, dt);
    
    % 7.apply boundary conditions to ghost cell velocities only using initial time
    % u = bcGhost_u(u, t);
    % v = bcGhost_v(v, t);

% 8. Start time loop-------------------------------------

ifig = 1; % figure counter
numTimeSteps = 0; %number of time steps counter

for k = 1:length(times_of_interest)
    
    %set target time
    endTime = times_of_interest(k); 

    while (t < endTime)
    %1. find dt using Adam bashforth limitation
    dt  = calcDtBurgers561(t, endTime, u, v);
        
    %2. use CN to solve Stokes Momentum eq w/o pressure
        % we will first use adam-bashforth to calculate the hyperbolic
        % terms for u and v, add those to the souce terms and then use CN
        % to find u and v. Then we will repeat this process, but for the Y
        % terms

        %record previous u and v terms
        u_prev = u;
        v_prev = v;

        %Find q at t_n
        [Qu, Qv, QY_n] = calcSourceIB(u,v,Y,t,dt);

        %Hyperbolic terms will be added to source terms to "freeze the hyperbolic term in time during ADI"
        [Hu,Hv] = hyperbolic_uv_2D(u,v);
        Qu = Qu + (1.5)*Hu - (0.5)*Hu_prev;
        Qv = Qv + (1.5)*Hv - (0.5)*Hv_prev;       
        
        %perform step one of ADI CN
        u = parabolic_CN1_2D_u(u, Qu, dt);
        v = parabolic_CN1_2D_v(v, Qv, dt);

        %Find q at t_n+1/2
        [Qu, Qv, QY_n12] = calcSourceIB(u,v,Y,t+dt/2,dt);
        Qu = Qu + (1.5)*Hu - (0.5)*Hu_prev;
        Qv = Qv + (1.5)*Hv - (0.5)*Hv_prev;

        %perform step 2 of ADI CN
        u = parabolic_CN2_2D_u(u, Qu, dt);
        v = parabolic_CN2_2D_v(v, Qv, dt);
        
        %%calculate Y parabolic terms    
        HY = hyperbolic_Y_WENO_2D(Y, u_prev, v_prev, u, v, dt);
        QY = QY_n + (1.5)*HY - (0.5)*HY_prev; 
               
        Y = parabolic_CN1_2D_Y3(Y, QY, dt);
   
        QY = QY_n12 + (1.5)*HY - (0.5)*HY_prev; 

        Y = parabolic_CN2_2D_Y3(Y, QY, dt);
        
        %increment t by half time step
        t = t+dt;

        %update previous hyperbolic terms
        Hu_prev = Hu;
        Hv_prev = Hv;
        HY_prev = HY;

   % 3. apply boundary conditions to boundaries and ghost cells using tn+1
        u = bc_u(u, t);
        v = bc_v(v, t);
        Y = bc_Y3(Y, t);

   % 4. correct outlet velocities to ensure volume conservation
        % [u, v] = correctOutlet(u,v);

   % 5. calculate right hand side of Poisson equation
        % f= (1/dt)*calcDivV(u, v);

   % 6. solve Poisson equation for Lagrange multiplier using Neumann boundary conditions
        % phi = bcGS(phi); %note we are using the same phi as in the last step... should reduce number of iters needed
        % [phi, Linf, iter] = myPoisson(phi, f, h, nIterMax, epsilon);

   % 7.project velocities using Lagrange multiplier
        % [u, v]= projectV(u, v, phi, dt);

   % 8. apply boundary conditions to ghost cell velocities only using tn+1
        % u = bcGhost_u(u, t);
        % v = bcGhost_v(v, t);

   %9. repeat till time step reached
    
    %increment number of time steps taken
    numTimeSteps = numTimeSteps+1;

   %calculate k and s
    k_values(numTimeSteps) = calck(u, v);
    S_values(numTimeSteps) = calcS3(Y);

    %store time
    Time_values(numTimeSteps) = t;
    
    %check if k, s and times are going to be too short --> if so increase
    %preallocated memory size
     if( numTimeSteps+(BlockSize/10) > currentSize )  % less than 10%*BLOCK_SIZE free slots
        currentSize = currentSize + BlockSize;       % add new BLOCK_SIZE slots
        k_values(numTimeSteps+1:currentSize) = 0;
        S_values(numTimeSteps+1:currentSize) = 0;
        Time_values(numTimeSteps+1:currentSize) = 0;
     end

    end

    if(ismember(endTime,Contours_times_of_interest ))
        %plot u
        examFig1 = figure(1);
        subplot(3,3,ifig);
        pcolor(xf, yc, u')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([-2.5, 2.5]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('u at time =', num2str(endTime)));
        colormap jet;
        colorbar;
        
         %plot v
        examFig2 = figure(2);
        subplot(3,3,ifig);
        pcolor(xc, yf, v')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([-2.5, 2.5]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('v at time =', num2str(endTime)));
        colormap jet;
        colorbar;
        
         %plot Y
        examFig3 = figure(3);
        subplot(3,3,ifig);
        pcolor(xc3, yc3, Y')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([0 1]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('Y at time =', num2str(endTime)));
        colormap hot;
        colorbar;

        %increment counter
        ifig = ifig+1;
    end
    
end

%plot k and S
examFig4 = figure(4);
plot(Time_values(1:numTimeSteps), k_values(1:numTimeSteps));
set(gca,'fontsize',16);
xlabel('time'); ylabel('k(t)');
title(' k(t) vs t');

examFig5 = figure(5);
plot(Time_values(1:numTimeSteps), S_values(1:numTimeSteps));
set(gca,'fontsize',16);
xlabel('time'); ylabel('S');
title('S vs t');
ylim([0, .25]);


K = 0.1 * myTrapezoidal(Time_values(1:numTimeSteps), k_values(1:numTimeSteps))
R = 0.1 * myTrapezoidal(Time_values(1:numTimeSteps), S_values(1:numTimeSteps))
