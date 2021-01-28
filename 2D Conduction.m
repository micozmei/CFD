% 2D Heat Equation
clear; close all; clc

% Parameters
alpha = 1;
L  = 1;
dx = 0.05;
dy = dx;
n  = L/dx+1;    
x  = 0:dx:L;
y  = x;
dt = 6.25e-4;
t  = dt*50;
nt = t/dt;

% Tolerance for Steady State
% TOL = 1e-6; 

% Initial Conditions
T = zeros(n);
T0 = zeros(1,n);
for i=1:n
    if x(i)<=0.5
        T0(i)=2*x(i);
    else
        T0(i)=2-2*x(i);
    end
end
T(1,1:n) = T0; % bottom
T(n,1:n) = 0;  % top
T(1:n,1) = 0;  % left
T(1:n,n) = 0;  % right

% Iterations for Steady State
% error = 1; 
% k = 0;

% while error > TOL
%   k = k+1;
for k = 1:nt 
    T_old = T;
    % Boundary Conditions are constant so skip indices 1 and n in the loop
    for i = 2:n-1
        for j = 2:n-1
            T(i,j) = alpha*dt*((T_old(i+1,j)-2*T_old(i,j)+T_old(i-1,j))/dx^2 ...
                + (T_old(i,j+1)-2*T_old(i,j)+T_old(i,j-1))/dy^2) ...
                + T_old(i,j);
        end
    end
%   error = max(max(abs(T_old-T)));
%   convergence(k,:)=error;
end

% Plots
figure(1)
pcolor(x,y,T),shading interp
title(['\Deltat = ',num2str(dt),', ','t = ',num2str(t)],'FontSize',12')
xlabel('x'),ylabel('y'),colorbar

% figure()
% plot(convergence)
% title('Convergence')
% xlabel('Iteration'),ylabel('Error'),xlim([0 100])
%
% figure()
% pcolor(x,y,T),shading interp
% title('Temperature (Steady State)')
% xlabel('x'),ylabel('y'),colorbar
% 
% figure()
% surf(x,y,T)
% title('Temperature (Steady State)')
% xlabel('x'),ylabel('y'),zlabel('T'),colorbar,view(3)