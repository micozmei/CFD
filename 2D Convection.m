% 2D Convection Equation
clear; close all; clc

% Parameters
alpha_x = 0.5;
alpha_y = 0.5;
L  = 20;
dx = 0.05;
dy = dx;
nx = 2*L/dx+1;
ny = 2*L/dy+1;
x  = -L:dx:L;
y  = x;
CFL = 0.6;
dt = CFL/(alpha_x/dx+alpha_y/dy);
t  = 15;
nt = t/dt;

% Initial Conditions
T = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        if ((-7+alpha_y*t<=y(j))&&(y(j)<=-5+alpha_y*t)&&...
            (-7+alpha_x*t<=x(i))&&(x(i)<=-5+alpha_x*t))
            T(j,i)=1;
        else
            T(j,i)=0;
        end
    end
end

% Boundary Conditions
T(1,:)  = 0; % bottom
T(nx,:) = 0; % top
T(:,1)  = 0; % left
T(:,ny) = 0; % right

% Numercial Scheme
for k = 1:nt 
    T_old = T;
    for i = 2:nx-1
        for j = 2:ny-1
            % Explicit Upwind
            T(i,j) = -alpha_x*dt*((T_old(i,j)-T_old(i-1,j))/dx) ...
                     -alpha_y*dt*((T_old(i,j)-T_old(i,j-1))/dy) ...
                     +T_old(i,j);
            % Lax-Wendroff
%             T(i,j) = -alpha_x*dt*((T_old(i+1,j)-T_old(i-1,j))/(2*dx)) ...
%                      -alpha_y*dt*((T_old(i,j+1)-T_old(i,j-1))/(2*dy)) ...
%                      +((alpha_x*dt)^2)*((T_old(i+1,j)-2*T_old(i,j)+T_old(i-1,j))/(2*dx*dx)) ...
%                      +((alpha_y*dt)^2)*((T_old(i,j+1)-2*T_old(i,j)+T_old(i,j-1))/(2*dy*dy)) ...
%                      +((dt^2)/(8*dx*dy))*(2*alpha_x*alpha_y)*((T_old(i+1,j+1)-T_old(i-1,j+1))-(T_old(i+1,j-1)-T_old(i-1,j-1)))...
%                      +T_old(i,j);
        end
    end
end

% Plots
figure(1)
pcolor(x,y,T),shading interp
title(['t = 15, ','CFL = ',num2str(CFL)],'FontSize',12')
xlabel('x'),ylabel('y'),colorbar

figure(2)
surf(x,y,T),shading interp
title(['t = 15, ','CFL = ',num2str(CFL)],'FontSize',12')
xlabel('x'),ylabel('y'),zlabel('T'),colorbar,view(3)