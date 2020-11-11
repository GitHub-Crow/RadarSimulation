FileName = 'D:/Assignment/Workspace/Radar Assignments/Materials/Chapter3Material/hw4.dat';
fid = fopen(FileName);
data = textscan(fid, '%f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
data = cell2mat(data);
permittivity = 80 - 17*1j;
lambda = 107;
k = 2*pi/lambda;

D = data(:, 1);
sa_back =  data(:, 2).^2 + data(:, 3).^2;
sb_back = data(:, 4).^2 + data(:, 5).^2;
sa_forward = data(:, 6).^2 + data(:, 7).^2;
sb_forward = data(:, 8).^2 + data(:, 9).^2;

aspect_r = 0.9951 + 0.0251.*D - 0.03644*D.^2 + 0.005303*D.^3 - 0.0002492*D.^4;
g_squared = 1./(aspect_r).^2 - 1;
Lb = (1 + g_squared)./g_squared.*(1 - 1./(sqrt(g_squared)).*atan(sqrt(g_squared)));
La = (1 - Lb)./2;

% Znric
% sa_re = abs(pi^2.*D.^3/(6*lambda^2).*(La + 1/(permittivity - 1)));
% sb_re = abs(pi^2.*D.^3/(6*lambda^2).*(Lb + 1/(permittivity) - 1));


% Guifu Zhang
sa_re = pi/2*k^2.*(D/2).^3.*(real(permittivity) - 1)./(3*(1 + (real(permittivity) - 1)*La));
sb_re = pi/2*k^2.*(D/2).^3.*(real(permittivity) - 1)./(3*(1 + (real(permittivity) - 1)*Lb));

figure(1);
semilogy(D, sa_back, 'Color', 'r',  'LineStyle', '-', 'LineWidth', 2);
hold on
semilogy(D, sa_forward, 'Color', 'r',  'LineStyle', '--', 'LineWidth', 2);
hold on
semilogy(D, sb_back, 'Color', 'r',  'LineStyle', '-.', 'LineWidth', 2);
hold on
semilogy(D, sb_forward, 'Color', 'r',  'LineStyle', ':', 'LineWidth', 2);
hold on
semilogy(D, sa_re, 'Color', 'k',  'LineStyle', '-', 'LineWidth', 2);
hold on
plot(D, sb_re, 'Color', 'k',  'LineStyle', ':', 'LineWidth', 2);


legend('real backscatter S_a', 'real forwardscatter S_a','real backscatter S_b','real forwardscatter S_a',...
       'Rayleigh scatter S_a','Rayleigh scatter S_b',...
       'Location','southeast');