%% MAE 561-- Final Project

clear
clear all

clf reset

global h CFL Re Sc xf xc yf yc yc3 xc3 Lx Ly t;

%Declare some constants
M=192;
N=128;
CFL = 0.9;
Re = 500;
Sc = 20;
Lx=3;
Ly=2;

times_of_interest = .1:.1:10;

%preallocate space for movie
loops = length(times_of_interest);
YMov(loops) = struct('cdata',[],'colormap',[]);
uMov(loops) = struct('cdata',[],'colormap',[]);
vMov(loops) = struct('cdata',[],'colormap',[]);
quivMov(loops) = struct('cdata',[],'colormap',[]);

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

% 3. correct outlet velocities to ensure volume conservation
[u, v] = correctOutlet(u,v);

% 4. calculate right hand side of Poisson equation using ∆t = 1
dt = 1;
f = (1/dt)*calcDivV(u, v);

% 5. solve Poisson equation for Lagrange multiplier using Neumann boundary conditions
phi = ones(M+2, N+2); %inialized phi = 1
phi = bcGS(phi);
[phi, Linf, iter] = myPoisson(phi, f, h, nIterMax, epsilon);

% 6. project velocities using Lagrange multiplier using ∆t = 1
[u, v]= projectV(u, v, phi, dt);

% 7.apply boundary conditions to ghost cell velocities only using initial time
u = bcGhost_u(u, t);
v = bcGhost_v(v, t);

%Adams-bashforth requires hyperbolic terms from previous step -->
%initialize these previous terms for n-1 using values at n
dt  = calcDtBurgers561(t, times_of_interest(2), u, v);
HY_prev = hyperbolic_Y_WENO_2D(Y, u, v, u, v, dt); 
[Hu_prev,Hv_prev] = hyperbolic_uv_2D(u,v);

% 8. Start time loop-------------------------------------

numTimeSteps =1; %number of time steps counter
outputFlag = 0;

%calculate initial k and S
%calculate k and s
k_values(1) = calck(u, v);
S_values(1) = calcS3(Y);
Time_values(1) = t;

for k = 1:length(times_of_interest)
    
    while(~outputFlag)
    %1. find dt using Adam bashforth limitation
    [dt, outputFlag]  = calcDtBurgers561(t, times_of_interest(k), u, v);
        
    %2. use CN to solve Stokes Momentum eq w/o pressure
        % we will first use adam-bashforth to calculate the hyperbolic
        % terms for u and v, add those to the souce terms and then use CN
        % to find u and v. Then we will repeat this process, but for the Y
        % terms

        %record previous u and v terms
        u_prev = u;
        v_prev = v;

        %Find q at t_n
        [Qu, Qv, QY_n] = calcSourceIBFinal(u,v,Y,t,dt);

        %Hyperbolic terms will be added to source terms to "freeze the hyperbolic term in time during ADI"
        [Hu,Hv] = hyperbolic_uv_2D(u,v);
        Qu = Qu + (1.5)*Hu - (0.5)*Hu_prev;
        Qv = Qv + (1.5)*Hv - (0.5)*Hv_prev;       
        
        %perform step one of ADI CN
        u = parabolic_CN1_2D_u(u, Qu, dt);
        v = parabolic_CN1_2D_v(v, Qv, dt);

        %Find q at t_n+1/2
        [Qu, Qv, QY_n12] = calcSourceIBFinal(u,v,Y,t+dt/2,dt);
        Qu = Qu + (1.5)*Hu - (0.5)*Hu_prev;
        Qv = Qv + (1.5)*Hv - (0.5)*Hv_prev;

        %perform step 2 of ADI CN
        u = parabolic_CN2_2D_u(u, Qu, dt);
        v = parabolic_CN2_2D_v(v, Qv, dt);
        
        %%calculate Y parabolic terms    
        HY = hyperbolic_Y_WENO_2D(Y, u_prev, v_prev, u, v, dt);
        QY = QY_n + HY; 
               
        Y = parabolic_CN1_2D_Y3(Y, QY, dt);
   
        QY = QY_n12 + HY; 

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
        [u, v] = correctOutlet(u,v);

   % 5. calculate right hand side of Poisson equation
        f= (1/dt)*calcDivV(u, v);

   % 6. solve Poisson equation for Lagrange multiplier using Neumann boundary conditions
        phi = bcGS(phi); %note we are using the same phi as in the last step... should reduce number of iters needed
        [phi, Linf, iter] = myPoisson(phi, f, h, nIterMax, epsilon);

   % 7.project velocities using Lagrange multiplier
        [u, v]= projectV(u, v, phi, dt);

   % 8. apply boundary conditions to ghost cell velocities only using tn+1
        u = bcGhost_u(u, t);
        v = bcGhost_v(v, t);

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
    
    outputFlag = 0;


    
        %plot u
        examFig1 = figure(1);
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        pcolor(xf, yc, u')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([-2.5, 2.5]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('u at time =', num2str(times_of_interest(k))));
        colormap jet;
        colorbar;
        drawnow
        uMov(k) = getframe(examFig1);

         %plot v
        examFig2 = figure(2);
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        pcolor(xc, yf, v')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([-2.5, 2.5]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('v at time =', num2str(times_of_interest(k))));
        colormap jet;
        colorbar;
        drawnow
        vMov(k) = getframe(examFig2);
        
        %plot Y
        examFig3 = figure(3);
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        pcolor(xc3, yc3, Y')
        shading interp;
        ylim([0,2]);
        xlim([0,3]);
        clim([0 1]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('Y at time =', num2str(times_of_interest(k))));
        colormap hot;
        colorbar;
        drawnow
        YMov(k) = getframe(examFig3);

        examFig4 = figure(4);
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        pcolor(xc3, yc3, Y')
        shading interp;
        clim([0 1]);
        colormap hot;
        colorbar;
        hold on
        quiv = quiver(xf,yf, u(1:M+1,1:N+1)',v(1:M+1,1:N+1)',2);
        % [X_mesh, Y_mesh] = meshgrid(xf,yf);
        % verts = stream2(X_mesh,Y_mesh,u(1:M+1,1:N+1)',v(1:M+1,1:N+1)',3*ones([length(yf), 1]), yf);
        % streamline(verts)
        ylim([0,2]);
        xlim([0,3]);
        set(gca,'fontsize',16);
        xlabel('x'); ylabel('y');
        title(strcat('Quivers at time =', num2str(times_of_interest(k))));
        hold off
        drawnow
        quivMov(k) = getframe(examFig4);

end

fig1 = figure
movie(fig1, uMov);

fig2 = figure
movie(fig2, vMov);

fig3 = figure
movie(fig3, YMov);

fig4 = figure
movie(fig4, quivMov);



% K = 0.1 * myTrapezoidal(Time_values(1:numTimeSteps), k_values(1:numTimeSteps))
% R = 0.1 * myTrapezoidal(Time_values(1:numTimeSteps), S_values(1:numTimeSteps))