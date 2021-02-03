%% Solution to Sod's Shock Tube Problem
clear all; close all; clc
global P_l rho_l P_r rho_r gamma mu Beta

% |-- L --|-- 1 --|-- 2 --|-- 3 --|-- R --|
% "post" refers to stations 2 and 3

% Left
% x1 = Start of the rarefaction wave that is moving to the left
% x2 = End of the rarefaction wave that is moving to the left;
% x0 = Location of the initial membrane
% x3 = Location of the contact discontinuity separating left/right fluid
% x4 = Location of the shock front that is moving to the right
% Right

% Initial Conditions
x0    = 0;
P_l   = 1;
rho_l = 1;
u_l   = 0;
s_l   = 0;
P_r   = 0.1;
rho_r = 0.125;
u_r   = 0;
gamma = 1.4;
mu    = sqrt((gamma-1)/(gamma+1));
Beta  = (gamma-1)/(2*gamma);
R     = 1; % Gas Constant

% Speed of Sound
c_l        = power((gamma*P_l/rho_l),0.5);
M_l        = u_l/c_l;
c_r        = power((gamma*P_r/rho_r),0.5);
M_r        = u_r/c_r;
P_post     = fzero(@iterative_calc,3); % Modify Initial Guess
v_post     = 2*(sqrt(gamma)/(gamma-1))*(1-power(P_post,(gamma-1)/(2*gamma)));
rho_post   = rho_r*(((P_post/P_r)+ mu^2)/(1 + mu*mu*(P_post/P_r)));
v_shock    = v_post*((rho_post/rho_r)/((rho_post/rho_r)-1));
rho_middle = (rho_l)*power((P_post/P_l),1/gamma);

% Boundaries (can be set)
x_min    = -0.5;
x_max    = 0.5;

% Time
t = 0.75*(x_max/v_shock); % Set to time for shock to reach 75% of x_max

% Key Positions
x1 = x0 - c_l*t;
x3 = x0 + v_post*t;
x4 = x0 + v_shock*t;

% Determining x2
c_2 = c_l - ((gamma-1)/2)*v_post;
M_2 = v_post/c_2;
x2  = x0 + (v_post-c_2)*t;

% Data Structure
n_points = 10000;
x_a      = linspace(x_min,x_max,n_points);
data.x   = x_a';
data.rho = zeros(n_points,1); % Density
data.P   = zeros(n_points,1); % Pressure
data.u   = zeros(n_points,1); % Velocity
data.c   = zeros(n_points,1); % Speed of Sound
data.M   = zeros(n_points,1); % Mach Number
data.e   = zeros(n_points,1); % Specific Internal Energy
data.E   = zeros(n_points,1); % Specific Total Energy
data.H   = zeros(n_points,1); % Specific Enthalpy
data.s   = zeros(n_points,1); % Specific Entropy

for index = 1:n_points
    if data.x(index) < x1
        % Solution before x1
        data.rho(index) = rho_l;
        data.P(index)   = P_l;
        data.u(index)   = u_l;
        data.c(index)   = c_l;
        data.M(index)   = M_l;
    elseif (x1 <= data.x(index) && data.x(index) <= x2)
        % Solution between x1 and x2
        data.c(index)   = mu*mu*((x0 - data.x(index))/t) + (1-mu*mu)*c_l;
        data.rho(index) = rho_l*power((data.c(index)/c_l),2/(gamma-1));
        data.P(index)   = P_l*power((data.rho(index)/rho_l),gamma);
        data.u(index)   = (1-mu*mu)*((-(x0-data.x(index))/t)+c_l);
        data.M(index)   = data.u(index)/data.c(index);
    elseif (x2 <= data.x(index) && data.x(index) <= x3)
        % Solution between x2 and x3
        data.rho(index) = rho_middle;
        data.P(index)   = P_post;
        data.u(index)   = v_post;
        data.c(index)   = c_2;
        data.M(index)   = M_2;
    elseif (x3 <= data.x(index) && data.x(index) <= x4)
        % Solution between x3 and x4
        data.rho(index) = rho_post;
        data.P(index)   = P_post;
        data.u(index)   = v_post;
        data.c(index)   = sqrt(gamma*data.P(index)/data.rho(index));
        data.M(index)   = data.u(index)/data.c(index);
        s_3             = s_l + R*log(((P_post/P_post)^(1/(gamma-1)))*...
                          ((rho_post/rho_middle)^(-gamma/(gamma-1))));
        data.s(index)   = s_3;
    elseif x4 < data.x(index)
        % Solution after x4
        data.rho(index) = rho_r;
        data.P(index)   = P_r;
        data.u(index)   = u_r;
        data.c(index)   = c_r;
        data.M(index)   = M_r;
        data.s(index)   = s_3 + R*log(((P_r/P_post)^(1/(gamma-1)))*...
                          ((rho_r/rho_post)^(-gamma/(gamma-1))));
    end
    data.e(index) = data.P(index)/((gamma-1)*data.rho(index));
    data.E(index) = data.e(index)+0.5*(data.u(index))^2;
    data.H(index) = data.E(index)+data.P(index)/data.rho(index);
end

%% MacCormack Solution to Sod's Shock Tube Problem
% Initialization
nx   = 201; % Number of grid points
rho  = zeros(nx,1);
u    = zeros(nx,1);
E    = zeros(nx,1);
P    = zeros(nx,1);
Q    = zeros(nx,3);
F    = zeros(nx,3);
Qnew = zeros(nx,3);
c    = zeros(nx,1);
s    = zeros(nx,1);

% Parameters
CFL   = 0.7;
dx    = 1/(nx-1);

% Create Mesh
x = (-0.5:dx:0.5);

% Initial Conditions
for i=1:nx
    if (x(i) > 0)
        P(i)   = 0.1;
        rho(i) = 0.125;
        u(i)   = 0.0;
    elseif (x(i) <= 0)
        P(i)   = 1.0;
        rho(i) = 1.0;
        u(i)   = 0.0;
    end
    E(i) = P(i)./(gamma-1) + 0.5*rho(i)*(u(i).^2);
    c(i) = sqrt(gamma*P(i)/rho(i));
end

% Qnew is the solution vector
Qnew(:,1) = rho;
Qnew(:,2) = rho.*u;
Qnew(:,3) = E;

t_current=0;
% Time Loop
while (t_current<t)
    dt = CFL*dx./(max(max(c+sqrt(u.*u))));
    for i=1:nx
        % Flux Calculation
        Q(i,1) = rho(i);
        Q(i,2) = rho(i).*u(i);
        Q(i,3) = E(i);

        F(i,1) = rho(i).*u(i);
        F(i,2) = rho(i).*u(i).^2 + P(i);
        F(i,3) = u(i).*(E(i) + P(i));
    end

    for i=2:nx-1
        % Predictor Step
        Qstar(i,1) = Q(i,1) - (dt/dx)*(F(i+1,1) - F(i,1));
        Qstar(i,2) = Q(i,2) - (dt/dx)*(F(i+1,2) - F(i,2));
        Qstar(i,3) = Q(i,3) - (dt/dx)*(F(i+1,3) - F(i,3));

        % Corrector Step
        ustar(i) = Qstar(i,2) / Qstar(i,1);
        Pstar(i) = (gamma-1)*(Qstar(i,3) - 0.5*(Qstar(i,2).^2/Qstar(i,1)));

        Fstar(i,1) = Qstar(i,1).*ustar(i);
        Fstar(i,2) = Qstar(i,1).*(ustar(i).^2) + Pstar(i);
        Fstar(i,3) = ustar(i).*(Qstar(i,3) + Pstar(i));
        Fstar(1,:) = Fstar(2,:); % Boundary Condition
        
        Qnew(i,1) = 0.5*(Qnew(i,1) + Qstar(i,1) - (dt/dx) * (Fstar(i,1) - Fstar(i-1,1)));
        Qnew(i,2) = 0.5*(Qnew(i,2) + Qstar(i,2) - (dt/dx) * (Fstar(i,2) - Fstar(i-1,2)));
        Qnew(i,3) = 0.5*(Qnew(i,3) + Qstar(i,3) - (dt/dx) * (Fstar(i,3) - Fstar(i-1,3)));
    end
    
    % Boundary Conditions
    Qnew(1,:)   = Qnew(2,:);
    Qnew(nx,:)  = Qnew(nx-1,:);
    Fstar(nx,:) = Fstar(nx-1,:);

    % Primitive Variables
    rho = Qnew(:,1);
    u   = Qnew(:,2)./rho;
    E   = Qnew(:,3);
    P   = (gamma-1)*(E-0.5*rho.*(u.*u));
    c   = sqrt(gamma*P./rho);
    
    % Advance in Time
    t_current = t_current+dt; 
end

M = u./c;
E = E./rho;
e = E-0.5*u.^2;
H = E+P./rho;

s(1) = s_l;
for i=2:length(x)
    s(i) = s(i-1) + R*log(((P(i)/P(i-1))^(1/(gamma-1)))*...
                          ((rho(i)/rho(i-1))^(-gamma/(gamma-1))));
end

%% Plots
figure()
plot(x_a,data.rho,'-k','LineWidth',2),hold on
plot(x,rho,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Density (kg/m^3)')
title('Density vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthEast')

figure()
plot(x_a,data.P,'-k','LineWidth',2),hold on
plot(x,P,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Pressure (Pa)')
title('Pressure vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthEast')

figure()
plot(x_a,data.u,'-k','LineWidth',2),hold on
plot(x,u,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Velocity (m/s)')
title('Velocity vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.c,'-k','LineWidth',2),hold on
plot(x,c,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Speed of Sound (m/s)')
title('Speed of Sound vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.M,'-k','LineWidth',2),hold on
plot(x,M,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Mach Number')
title('Mach Number vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.e,'-k','LineWidth',2),hold on
plot(x,e,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Specific Internal Energy (J/kg)')
title('Specific Internal Energy vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.E,'-k','LineWidth',2),hold on
plot(x,E,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Specific Total Energy (J/kg)')
title('Specific Total Energy vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.H,'-k','LineWidth',2),hold on
plot(x,H,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Specific Enthalpy (J/kg)')
title('Specific Enthalpy vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

figure()
plot(x_a,data.s,'-k','LineWidth',2),hold on
plot(x,s,'or','LineWidth',1)
xlabel('x (m)'),ylabel('Specific Entropy (J/kg K)')
title('Specific Entropy Change vs Position'),grid on
legend('Analytical Solution','MacCormack Solution','Location','NorthWest')

% Iteratively solve for Pressure
function y = iterative_calc(P)
global P_l rho_l P_r rho_r gamma mu Beta
y = (P_l^(Beta)-P^Beta)*sqrt((1-mu^4)*(P_l^(1/gamma))/(rho_l*mu^4))-...
    (P-P_r)*sqrt((1-mu^2)/(rho_r*(P+P_r*mu^2)));
end