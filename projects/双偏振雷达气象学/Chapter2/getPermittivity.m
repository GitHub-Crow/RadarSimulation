function [constant] = getPermittivity(f, epsilon_s, epsilon_inf, alpha, lambda_s, delta)
c = 3*10^8;
epsilon0 = 8.854*10^-12;
lambda = c/f;
real = epsilon_inf + (epsilon_s - epsilon_inf)*(1 + (lambda_s/lambda)^(1-alpha)*sin(alpha*pi/2))/(1 + 2*(lambda_s/lambda)^(1-alpha)*sin(alpha*pi/2) + (lambda_s/lambda)^(2-2*alpha));
imag = (epsilon_s - epsilon_inf)*(lambda_s/lambda)^(1-alpha)*cos(alpha*pi/2)/(1 + 2*(lambda_s/lambda)^(1-alpha)*sin(alpha*pi/2) + (lambda_s/lambda)^(2-2*alpha)) + delta*lambda/(2*pi*c*epsilon0);
constant = real - imag*1j;
end

