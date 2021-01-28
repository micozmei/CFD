% 1D Convection Equation
clear all; close all; clc

alpha=0.5; % Convection Coefficient
L=20;      % Length of Rod

dx=0.05;       % Position Step
x=-L/2:dx:L/2; % Position Intervals
N=L/dx;        % Number of Points in Space

CFL=0.6;         % CFL = alpha*dt/dx
dt=CFL*dx/alpha; % Time Step

T=15;   % Time Interval
M=T/dt; % Number of Time Steps

%% Analytical Solution of 1D Convection Equation 
u = zeros(length(x),1);
for i = 1:length(x)
    if x(i) < -7+alpha*T
        u(i)=0.0;
    elseif x(i) <= -5+alpha*T
        u(i)=1.0;
    else
        u(i)=0.0;
    end
end

%% Forward First-Order Upwind Solution of 1D Convection Equation 
% Initial Conditions (t=0)
u0_up = zeros(length(x),1); % Pre-allocate
for i=1:N+1
    if x(i) < -7
        u0_up(i)=0;
    elseif x(i) <= -5
        u0_up(i)=1;
    else
        u0_up(i)=0;
    end
end

% Partial Difference Equation (Numerical Scheme)
u1_up = u0_up; % Pre-allocate
u1_up_diff = u0_up;
u1_up_disp = u0_up;

for j=1:M % Time Iterations
    for i=2:N % Position Iterations
        u1_up(i)=u0_up(i)-CFL*(u0_up(i)-u0_up(i-1));
        % Diffusion Error (use central difference for d^2u/dx^2)
        d2u_dx2 = (u0_up(i+1)-2*u0_up(i)+u0_up(i-1))/(dx^2);
        u1_up_diff(i) = (alpha*dx/2)*(1-CFL)*d2u_dx2;
        % Dispersion Error (use central difference for d^3u/dx^3)
        if i > 2 && i < N
            d3u_dx3 = (u0_up(i+2)-2*u0_up(i+1)+2*u0_up(i-1)-u0_up(i-2))/(2*dx^3);
            u1_up_disp(i) = (alpha*dx*dx/6)*(3*CFL-2*CFL^2-1)*d3u_dx3;
        end
    end
    % Boundary Conditions (x=0)
    u1_up(1)=0;
    u0_up=u1_up; % Update BC
end

%% Lax-Wendroff Solution of 1D Convection Equation 
% Initial Conditions (t=0)
u0_lw = zeros(length(x),1); % Pre-allocate
for i=1:N+1
    if x(i) < -7
        u0_lw(i)=0;
    elseif x(i) <= -5
        u0_lw(i)=1;
    else
        u0_lw(i)=0;
    end
end

% Partial Difference Equation (Numerical Scheme)
u1_lw = u0_lw; % Pre-allocate
u1_lw_diff = u0_lw;
u1_lw_disp = u0_lw;

for j=1:M % Time Iterations
    for i=2:N % Position Iterations
        u1_lw(i)=u0_lw(i)-0.5*CFL*(u0_lw(i+1)-u0_lw(i-1))+0.5*CFL*CFL*...
                 (u0_lw(i+1)-2*u0_lw(i)+u0_lw(i-1));
        % Diffusion Error (use central difference for d^2u/dx^2)
        d2u_dx2 = (u0_lw(i+1)-2*u0_lw(i)+u0_lw(i-1))/(dx^2);
        u1_lw_diff(i) = 0;
        % Dispersion Error (use central difference for d^3u/dx^3)
        if i > 2 && i < N 
            d3u_dx3 = (u0_lw(i+2)-2*u0_lw(i+1)+2*u0_lw(i-1)-u0_lw(i-2))/(2*dx^3);
            u1_lw_disp(i) = (alpha*dx*dx/6)*(CFL^2-1)*d3u_dx3;
        end
    end
    % Boundary Conditions (x=0)
    u1_lw(1)=0;
    u0_lw=u1_lw; % Update BC
end

%% Implicit Euler Solution of 1D Convection Equation 
% Coefficients of the tridiagonal system
b = alpha/(2*dx); % Super diagonal: coefficients of u(i+1)
c = -b;           % Subdiagonal: coefficients of u(i-1)
a = 1/dt;         % Main Diagonal: coefficients of u(i)

% Tridiagonal matrix
A = diag(a*ones(1,N+1))+diag(b*ones(1,N),1)+diag(c*ones(1,N),-1);
% Fix coefficients of boundary nodes
A(1,1)=1; A(1,2)=0; A(end,end)=1; A(end,end-1)=0;

% Initial Conditions
u1_ie = zeros(length(x),1); % Pre-allocate
u1_ie_diff = u1_ie;
u1_ie_disp = u1_ie;
for i=1:N+1
    if x(i) < -7
        u1_ie(i)=0;
    elseif x(i) <= -5
        u1_ie(i)=1;
    else
        u1_ie(i)=0;
    end
end
% Boundary Conditions
ub=[0, 0];

% Loop over time steps
for j = 1:M
    d = [ub(1); u1_ie(2:N)/dt; ub(2)]; % Update LHS and preserve BCs
    u1_ie = A\d; % Solve the system
    for i = 2:N
        d2u_dx2 = (u1_ie(i+1)-2*u1_ie(i)+u1_ie(i-1))/(dx^2);
        u1_ie_diff(i) = (alpha*dx*CFL/2)*d2u_dx2;
        if i > 2 && i < N 
            d3u_dx3 = (u1_ie(i+2)-2*u1_ie(i+1)+2*u1_ie(i-1)-u1_ie(i-2))/(2*dx^3);
            u1_ie_disp(i) = -(alpha*dx*dx/6+alpha*dx*dx*CFL*CFL/3)*d3u_dx3;
        end
    end
end

%% Plot Solutions
figure(1)
plot(x,u1_up,x,u1_lw,x,u1_ie,x,u,'LineWidth',1.5)
title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
xlabel('x','FontSize',14),xlim([-10,10])
ylabel('u','FontSize',14),ylim([-0.4,1.4])
legend('Explicit Upwind','Lax-Wendroff','Implicit Euler',...
       'Analytical Solution','Location','NorthEast')

figure(2)
plot(x,u1_up_diff,x,u1_lw_diff,x,u1_ie_diff,'LineWidth',1.5)
title(['Diffusion Error (t = ',num2str(T),', ','CFL = ',num2str(CFL),')'],'FontSize',12)
xlabel('x','FontSize',14),xlim([-10,10])
ylabel('u','FontSize',14)
legend('Explicit Upwind','Lax-Wendroff',...
       'Implicit Euler','Location','NorthEast')
   
figure(3)
plot(x,u1_up_disp,x,u1_lw_disp,x,u1_ie_disp,'LineWidth',1.5)
title(['Dispersion Error (t = ',num2str(T),', ','CFL = ',num2str(CFL),')'],'FontSize',12)
xlabel('x','FontSize',14),xlim([-10,10])
ylabel('u','FontSize',14)
legend('Explicit Upwind','Lax-Wendroff',...
       'Implicit Euler','Location','NorthEast')
   
% figure(4)
% plot(x,u1_up_diff,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14)
% legend('Explicit Upwind Diffusion','Location','NorthWest')
% 
% figure(5)
% plot(x,u1_up_disp,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14)
% legend('Explicit Upwind Dispersion','Location','NorthWest')
% 
% figure(6)
% plot(x,u1_lw_diff,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14)
% legend('Lax-Wendroff Diffusion','Location','NorthWest')
% 
% figure(7)
% plot(x,u1_lw_disp,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14)
% legend('Lax-Wendroff Dispersion','Location','NorthWest')
% 
% figure(8)
% plot(x,u1_ie_diff,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14)
% legend('Implicit Euler Diffusion','Location','NorthWest')
% 
% figure(9)
% plot(x,u1_ie_disp,'LineWidth',1.5)
% title(['t = ',num2str(T),', ','CFL = ',num2str(CFL)],'FontSize',12)
% xlabel('x','FontSize',14),xlim([-10,10])
% ylabel('u','FontSize',14),ylim([-0.4,1.4])
% legend('Implicit Euler Dispersion','Location','NorthWest')