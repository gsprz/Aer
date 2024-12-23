clc


U_tilde = U_Inf_Mag;

span = config.Span;
CL_tilde =C_L; % @ 1.5°
AR = config.AspectRatio;

B1 = -CL_tilde/pi/config.AspectRatio;
theta=linspace(0,pi,200);
Gamma_theta = @(theta) 2*U_tilde*B1*span*sin(theta);
spanwise = @(theta) -span/2*cos(theta);

spanwise_norm = @(theta) spanwise(theta)./span;

CDi = pi*AR*B1^2

figure
plot(spanwise_norm(theta), abs(Gamma_theta(theta)))
grid on;
xlabel('apertura alare normalizzata');
ylabel('gamma');
title('Distribuzione di circolazione dell''ala ottima '); 

