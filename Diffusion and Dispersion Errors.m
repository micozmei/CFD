clear all; close all; clc
% Diffusion and Dispersion Errors
phi = linspace(0,pi,100);
CFL1 = 0.25;
CFL2 = 0.5;
CFL3 = 0.75;

eD_up_1 = sqrt(1-4.*CFL1.*(1-CFL1).*(sin(phi./2)).^2);
eD_up_2 = sqrt(1-4.*CFL2.*(1-CFL2).*(sin(phi./2)).^2);
eD_up_3 = sqrt(1-4.*CFL3.*(1-CFL3).*(sin(phi./2)).^2);
eP_up_1 = mod(atan(CFL1.*sin(phi)./(1-CFL1+CFL1.*cos(phi))),pi)./(CFL1.*phi);
eP_up_2 = mod(atan(CFL2.*sin(phi)./(1-CFL2+CFL2.*cos(phi))),pi)./(CFL2.*phi);
eP_up_3 = mod(atan(CFL3.*sin(phi)./(1-CFL3+CFL3.*cos(phi))),pi)./(CFL3.*phi);

eD_lw_1 = sqrt(1-4.*CFL1.*CFL1.*(1-CFL1.*CFL1).*((sin(phi./2)).^4));
eD_lw_2 = sqrt(1-4.*CFL2.*CFL2.*(1-CFL2.*CFL2).*((sin(phi./2)).^4));
eD_lw_3 = sqrt(1-4.*CFL3.*CFL3.*(1-CFL3.*CFL3).*((sin(phi./2)).^4));
eP_lw_1 = mod(atan(CFL1.*sin(phi)./(1-2.*CFL1.*CFL1.*(sin(phi/2)).^2)),pi)./(CFL1.*phi);
eP_lw_2 = mod(atan(CFL2.*sin(phi)./(1-2.*CFL2.*CFL2.*(sin(phi/2)).^2)),pi)./(CFL2.*phi);
eP_lw_3 = mod(atan(CFL3.*sin(phi)./(1-2.*CFL3.*CFL3.*(sin(phi/2)).^2)),pi)./(CFL3.*phi);

eD_ie_1 = sqrt(1./(1+CFL1.*CFL1.*(sin(phi)).^2));
eD_ie_2 = sqrt(1./(1+CFL2.*CFL2.*(sin(phi)).^2));
eD_ie_3 = sqrt(1./(1+CFL3.*CFL3.*(sin(phi)).^2));
eP_ie_1 = mod(atan(CFL1.*sin(phi)),pi)./(CFL1.*phi);
eP_ie_2 = mod(atan(CFL2.*sin(phi)),pi)./(CFL2.*phi);
eP_ie_3 = mod(atan(CFL3.*sin(phi)),pi)./(CFL3.*phi);

figure()
plot(phi,eD_up_1,phi,eD_up_2,phi,eD_up_3)
title('Explicit Upwind Diffusion Error')
xlabel('Phase Angle (rad)'),ylabel('Diffusion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])

figure()
plot(phi,eP_up_1,phi,eP_up_2,phi,eP_up_3)
title('Explicit Upwind Dispersion Error')
xlabel('Phase Angle (rad)'),ylabel('Dispersion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])

figure()
plot(phi,eD_lw_1,phi,eD_lw_2,phi,eD_lw_3)
title('Lax-Wendroff Diffusion Error')
xlabel('Phase Angle (rad)'),ylabel('Diffusion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])

figure()
plot(phi,eP_lw_1,phi,eP_lw_2,phi,eP_lw_3)
title('Lax-Wendroff Dispersion Error')
xlabel('Phase Angle (rad)'),ylabel('Dispersion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])

figure()
plot(phi,eD_ie_1,phi,eD_ie_2,phi,eD_ie_3)
title('Implicit Euler Diffusion Error')
xlabel('Phase Angle (rad)'),ylabel('Diffusion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])

figure()
plot(phi,eP_ie_1,phi,eP_ie_2,phi,eP_ie_3)
title('Implicit Euler Dispersion Error')
xlabel('Phase Angle (rad)'),ylabel('Dispersion Error')
legend('CFL = 0.25','CFL = 0.5','CFL = 0.75','Location','Best'),xlim([0,pi])