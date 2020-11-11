F = 2.8*10^9;
T = 0;
% 水的介电常数计算参数
w_epsilon_s = @(T) 78.54*(1 - 4.579*10^-3.*(T - 25) + 1.19*10^-5.*(T - 25)^2 - 2.8*10^-8.*(T - 25)^3);
w_epsilon_inf = @(T) 5.27137 + 2.16474*10^-2.*T - 1.31198*10^-3.*T.^2;
w_alpha = @(T) -16.8129./(T + 273) + 6.09265*10^-2;
w_lambda_s = @(T) 3.3836*10^-6*exp(2513.98./(T + 273));
w_delta = 1.1117*10^-4;
 disp(['water of permittivity at  ', num2str(T), ' : ', num2str(getPermittivity(F, w_epsilon_s(T), w_epsilon_inf(T), w_alpha(T), w_lambda_s(T), w_delta))]);


% 冰的介电常数计算参数
i_epsilon_s = @(T) 203.168 + 2.5.*T + 0.15.*T.^2;
i_epsilon_inf = @(T) 3.168;
i_alpha = @(T) 0.288 + 5.2*10^-3.*T + 2.3*10^-4.*T.^2;
i_lambda_s = @(T) 9.990288*10^-6*exp(6643.5./(T + 273));
i_delta = @(T) 1.1156*10^-13*exp(-6291.2./(T + 273));
permittivity_ice_0 = getPermittivity(F, i_epsilon_s(T), i_epsilon_inf(T), i_alpha(T), i_lambda_s(T), i_delta(T));
disp(['ice of permittivity at  ', num2str(T), ' : ', num2str(permittivity_ice_0)]);

% 利用M-G混合方程计算0度下雪的介电常数
length = 100;
density_snow = linspace(100, 917, length);
density_ice = 920;
permittivity_air = 1.0006;
y_air = (permittivity_ice_0 - permittivity_air)/(permittivity_ice_0 + 2*permittivity_air);
permittivity_snow_mg_air = @(ds) permittivity_air*(1 + 2.*(ds/density_ice)*y_air)./(1 - (ds/density_ice)*y_air);
y_ice = (permittivity_air - permittivity_ice_0)/(permittivity_air + 2*permittivity_ice_0);
permittivity_snow_mg_ice = @(ds) permittivity_ice_0*(1 + 2.*(1 - ds/density_ice)*y_ice)./(1 - (1 - ds/density_ice)*y_ice);
figure(1);
ps_mg_air = permittivity_snow_mg_air(density_snow);
subplot(2, 2, 1)
plot(density_snow, real(ps_mg_air));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Real of permittivity of snow given by M-G (air as background)');
subplot(2, 2, 2)
plot(density_snow, -imag(ps_mg_air));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Imag of permittivity of snow given by M-G (air as background)');
ps_mg_ice = permittivity_snow_mg_ice(density_snow);
subplot(2, 2, 3)
plot(density_snow, real(ps_mg_ice));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Real of permittivity of snow given by M-G (ice as background)');
subplot(2, 2, 4)
plot(density_snow, -imag(ps_mg_ice));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Imag of permittivity of snow given by M-G (ice as background)');
syms f p
eqn =  @(f, p) (1 - f).*(permittivity_air - p)./(permittivity_air + 2*p) + f.*(permittivity_ice_0 - p)./(permittivity_ice_0 + 2*p) == 0;
solution = solve(eqn, p);
ps_ps = subs(solution(2), f, density_snow/density_ice);
figure(2)
subplot(2, 1, 1)
plot(density_snow, real(ps_ps));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Real of permittivity of snow given by P-S');
subplot(2, 1, 2)
plot(density_snow, -imag(ps_ps));
xlabel('density of snow  /kg\cdotm^{-3}');
title('Imag of permittivity of snow given by P-S');

