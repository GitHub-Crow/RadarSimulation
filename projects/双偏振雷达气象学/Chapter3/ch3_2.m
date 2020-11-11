FileName = 'C:/space/Workspace/Radar/Materials/Chapter3Material/hm2p2.dat';
fid = fopen(FileName);
data = textscan(fid, '%f %f %f %f', 'HeaderLines', 1);
data = cell2mat(data);
permittivity = 41 - 41*1j;
lambda = 30;
k = 2*pi/lambda;
D = data(:, 1);
s_e = data(:, 2);
s_s = data(:, 3);
s_bs = data(:, 4);
s_g = pi/4.*D.^2;

% 瑞利散射效率因子
s_re_s = 8*pi*k^4/3*((real(permittivity) - 1) / (real(permittivity) + 2))^2.*(D./2).^6;
s_re_a = 4*pi*k/3*(-imag(permittivity))*(3 / (real(permittivity) + 2))^2.*(D./2).^3;
s_re_e = s_re_a + s_re_s;

% 米散射效率因子
[s_mie_e, s_mie_s] = getScatterAndExticationCoffByMie(D, lambda, permittivity);
s_mie_a = s_mie_e - s_mie_s;
figure(1)

plot(D, s_s./s_g, 'Color', 'r',  'LineStyle', '-', 'LineWidth', 2);
hold on 
plot(D, s_e./s_g, 'Color', 'r',  'LineStyle', '--', 'LineWidth', 2);
hold on 
plot(D, s_re_s./s_g, 'Color', 'g',  'LineStyle', '-', 'LineWidth', 2);
hold on 
plot(D, s_re_e./s_g, 'Color', 'g',  'LineStyle', '--', 'LineWidth', 2);
hold on 
plot(D, s_mie_s, 'Color', 'b',  'LineStyle', '-', 'LineWidth', 2);
hold on 
plot(D, s_mie_e, 'Color', 'b',  'LineStyle', '--', 'LineWidth', 2);
hold on 
axis([0, 200, 0.001, 10]);
legend('Real norm-Q_{s}', 'Real norm-Q_{t}',...
        'Rayleigh norm-Q_{s}', 'Rayleigh norm-Q_{t}',...
        'Mie norm-Q_{s}', 'Mie norm-Q_{t}');