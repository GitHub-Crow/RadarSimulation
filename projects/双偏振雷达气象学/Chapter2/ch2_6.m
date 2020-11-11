F = 2.8*10^9;
T = 0;

% 利用M-G混合方程计算0度下雪的介电常数
length = 100;
rw = linspace(0, 0.99, length);
density_drysnow = 100;
density_water = 1000;
density_snow = density_drysnow*(1 - rw.^2) + density_water*rw.^2;
density_ice = 920;
permittivity_air = 1.0006;
permittivity_water = 81 - 23.3*1j;
permittivity_ice = 3.17;
fi = (1 - rw).*density_snow/density_water;
fw = rw.*density_snow/density_water;
fa = 1 - fi - fw;

f = fa./(fa + fi);
y_ice = (permittivity_air - permittivity_ice)/(permittivity_air + 2*permittivity_ice);
ps_mg_hybird = permittivity_ice*(1 + 2.*f*y_ice)./(1 - f*y_ice);

f = 1 - fw;
y_water = (ps_mg_hybird - permittivity_water)./(ps_mg_hybird + 2*permittivity_water);
ps_mg_wetsnow = permittivity_water.*(1 + 2.*f.*y_water)./(1 - f.*y_water);

figure(1);
subplot(2, 1, 1)
plot(rw, real(ps_mg_wetsnow));
xlabel('ratio of melting');
title('Real of permittivity of snow given by M-G (water as background)');
subplot(2, 1, 2)
plot(rw, -imag(ps_mg_wetsnow));
xlabel('ratio of melting');
title('Imag of permittivity of snow given by M-G (water as background)');


