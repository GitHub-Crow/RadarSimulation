function calculateNC_W_R_Z(timeIndex, D, ND, V, ModelName)
NC = sum(ND, 1);
LWC = sum((pi/6*10^-6).*D.^3.*ND, 1);
R = sum(6*pi*10^-4.*D.^3.*V.*ND, 1);
Z = sum(D.^6.*ND, 1);  
figure();
subplot(4,1,1)
plot(timeIndex(:,1), NC);
title(strcat(ModelName, ':  Number Concentration  /m^{-3}'));
xlim([7 13])
subplot(4,1,2)
plot(timeIndex(:,1), LWC);
title(strcat(ModelName, ':  Liquid Water Content  /g\cdotm^{-3}'));
xlim([7 13])
subplot(4,1,3)
plot(timeIndex(:,1), R);
title(strcat(ModelName, ':  Rain Rate  /mm\cdoth^{-1}'));
xlim([7 13])
subplot(4,1,4)
plot(timeIndex(:,1), Z);
title(strcat(ModelName, ':  Reflectivity  /mm^{6}\cdotm^{-3}\cdotmm^{-1}'));
xlim([7 13])
end

